interface
SUBROUTINE COS_SZA(KSTART,KEND,KCOL,PGEMU,PGELAM,LDRADIATIONTIMESTEP,PMU0)
USE PARKIND1 , ONLY : JPIM, JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KSTART
INTEGER(KIND=JPIM),INTENT(IN) :: KEND
INTEGER(KIND=JPIM),INTENT(IN) :: KCOL
REAL(KIND=JPRB), INTENT(IN) :: PGEMU(KCOL)
REAL(KIND=JPRB), INTENT(IN) :: PGELAM(KCOL)
LOGICAL, INTENT(IN) :: LDRADIATIONTIMESTEP
REAL(KIND=JPRB), INTENT(OUT) :: PMU0(KCOL)
END SUBROUTINE COS_SZA
end interface