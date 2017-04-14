module subs
contains
subroutine check(status)
    USE NETCDF
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

end subroutine check

SUBROUTINE BilinealInterp(xP,yP,NETCDFValues,meltInterp, dxarg, dyarg, xInitarg, yInitarg)

        USE types
        USE CoordinateSystems
        USE SolverUtils
        USE ElementDescription
        USE DefUtils

        REAL(KIND=dp), INTENT(IN):: xP,yP
        REAL(Kind=dp), INTENT(IN) :: NETCDFValues(:,:)
        REAL(KIND=dp), INTENT(OUT) :: meltInterp
        REAL(KIND=dp) :: iFract, jFract, melt11, melt12, melt21, melt22
        INTEGER :: iIndex, jIndex
        REAL(KIND=dp) :: DX , DY, xInit, yInit
        REAL(KIND=dp), INTENT(IN), OPTIONAL ::  DXarg , DYarg, xInitarg, yInitarg


        if( .not. present(DXarg))then
            DX=2000.0
        else
            DX = DXarg
        end if

        if( .not. present(DYarg))then
            DY=2000.0
        else
            DY = DYarg
        end if

        if( .not. present(xInitarg))then
            xInit = -319000.0
        else 
            xInit = xInitarg
        end if

        if( .not. present(yInitarg))then
            yInit = -1000.0
        else
            yInit = yInitarg
        end if

        !if(associated(NETCDFValues)) then
        !        print *, 'associated'
        !else
        !        print *, 'not associated'
        !end if


        iIndex = int(floor((xP - xInit ) / DX)+1)
        iFract = (xP - xInit) - (iIndex-1) * DX

        jIndex = int(floor((yP - (yInit)) / DY)+1)
        jFract = (yP - (yInit)) - (jIndex-1) * DY

        if (iIndex .gt. 0) then
                melt11 = NETCDFValues(iIndex,jIndex)
                melt12 = NETCDFValues(iIndex,jIndex+1)
                melt21 = NETCDFValues(iIndex+1,jIndex)
                melt22 = NETCDFValues(iIndex+1,jIndex+1)
                meltInterp = 1/(DX*DY) * ( &
                        melt11*(DX-iFract)*(DY-jFract) + &
                        melt21*(iFract)*(DY-jFract) + &
                        melt12*(DX-iFract)*(jFract) + &
                        melt22*(iFract)*(jFract) &
                        )
        else    
                meltInterp = 0.0_dp
        end if

        if(xp.eq.458000.0_dp .and. yp.eq.41000.0_dp)then
                Print*, 'datas', iIndex, jIndex, iFract, jFract, melt11, melt12, melt21, melt22, meltInterp, NETCDFValues(iIndex,jIndex), NETCDFValues(iIndex,jIndex+1), NETCDFValues(iIndex+1,jIndex), NETCDFValues(iIndex+1,jIndex+1)
        end if

        !Print *, 'My PE: ', ParEnv % MyPE

        if (ParEnv % MyPE .eq. 14) then
                !Print *, NETCDFValues(iIndex,:)
        end if

End subroutine BilinealInterp
end module subs

SUBROUTINE MISOMIP_Melt( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils
  USE subs
  USE NETCDF

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Solver_t),POINTER :: PSolver
  TYPE(Variable_t),POINTER :: MeltVar=>NULL(), GMVar=>NULL()
 
  INTEGER , POINTER :: MeltPerm(:), GMPerm(:)
  REAL(KIND=dp) , POINTER :: Melt(:),GM(:)
  REAL(KIND=dp) , POINTER ::DATAPointer(:,:)
  INTEGER :: NMax, ncid, node, ncidDraft
  REAL(KIND=dp) ::  xP , yP, meltInt

  LOGICAL,SAVE :: Initialized = .FALSE.
  LOGICAL,SAVE :: ExtrudedMesh=.False.
  LOGICAL :: Found, Got

  CHARACTER(len = 200) :: FILE_NAME = '/scratch/cnt0021/gge6066/imerino/NEMO/WORK/nemo_ISOMIP/nemo_ISOMIP_EXP3_TEST/OUTPUT_1/ISOMIP-EXP3_TEST_1d_00000101_00000630_SBC.nc'
  CHARACTER(len = 200) :: FILE_NAME_DRAFT = '/scratch/cnt0021/gge6066/imerino/NEMO/WORK/nemo_ISOMIP/nemo_ISOMIP_EXP3_TEST/isf_draft_meter.nc'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='InitMELTMISOMIP'

  CHARACTER(len = 200), parameter :: meltname='fwfisf', xname='x'
  REAL(KIND=dp), allocatable, target :: meltvarNC(:,:), xVarNC(:)
  INTEGER :: varXid, varid, dimid1, dimid2, lenX, lenY, status1, res

!------------------------------------------------------------------------------
  Mesh => Model % Mesh

  NMax = Solver % Mesh % NumberOfNodes

!!! get required variables Zb,Zs,H

  MeltVar => VariableGet( Model % Mesh % Variables, 'Melt')
  IF (.NOT.ASSOCIATED(MeltVar)) THEN
     Message='Melt not found'
     CALL FATAL(SolverName,Message)
  END IF
  GMVar => VariableGet( Model % Mesh % Variables, 'GroundedMask')
  IF (.NOT.ASSOCIATED(GMVar)) THEN
     Message='GroundedMask not found'
     CALL FATAL(SolverName,Message)
  END IF

  FILE_NAME_DRAFT = GetString(Solver % Values,'Draft file',Got)
  FILE_NAME = GetString(Solver % Values,'Melt rates file',Got)

  Print *, 'filename', FILE_NAME_DRAFT

  MeltPerm => MeltVar % Perm
  Melt => MeltVar % Values

  GMPerm => GMVar % Perm
  GM => GMVar % Values

  CALL check(nf90_open(FILE_NAME,NF90_NOWRITE,ncid))
  CALL check(nf90_open(FILE_NAME_DRAFT,NF90_NOWRITE,ncidDraft))
  status1=nf90_inq_dimid(ncid,"x",dimid1)
  status1=nf90_inquire_dimension(ncid,dimid1,len=lenX)

  status1=nf90_inq_dimid(ncid,"y",dimid2)
  status1=nf90_inquire_dimension(ncid,dimid2,len=lenY)

  allocate(meltvarNC(lenX,lenY))
  allocate(DATAPOINTER(lenX,lenY))
  allocate(xVarNC(lenX))

  status1=nf90_inq_varid(ncid,meltname,varid)

  status1=nf90_get_var(ncid,varid,meltVarNC)

  status1=nf90_inq_varid(ncidDraft,xname,varXid)

  status1=nf90_get_var(ncidDraft,varXid,xVarNC)

  Print *, 'status ', status1

  DATAPOINTER => meltvarNC

  if (ParEnv % MyPE .eq. 14) then
         Print *, ' x=69 ', 'x = ', xVarNC(69), meltVarNC(69,:)
         Print *, ' x=69 ', 'x = ', xVarNC(70), meltVarNC(70,:)
         Print *, ' x=71 ', 'x = ', xVarNC(71), meltVarNC(71,:)
         Print *, ' x=72 ', meltVarNC(72,:)
         Print *, ' x=73 ', meltVarNC(73,:)
         Print *, ' x=74 ', meltVarNC(74,:)

  end if

  DO node=1, nMax 
        xP =  Mesh % Nodes % x(node)
        yP =  Mesh % Nodes % y(node)
        if (GM(GMPerm(node)) .lt. 0.0) then
                !Print *, 'hol', xP, yP, meltvarNC(80,20), DATAPointer(80,20), BiLinealInterp(xP,yP,DATAPointer)
                !Print *, 'hol', xP, yP, meltvarNC(80,20), DATAPointer(80,20)
                CALL BiLinealInterp(xP,yP,meltvarNC,meltInT)
                Melt(MeltPerm(node)) = meltInt * 1e-3 * 3600 * 24 * 365
                !Melt(MeltPerm(node)) = -100.0_dp
        else
                Melt(MeltPerm(node)) = 0.0_dp
        end if
  end do

  CALL check(nf90_close(ncid))
  CALL check(nf90_close(ncidDraft))

END SUBROUTINE MISOMIP_Melt

