SUBROUTINE rect_km(km,coord,e,v)
!
! This subroutine forms the "analytical" stiffness matrix for
! rectangular 4- or 8-node plane strain elements using nip=4.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::coord(:,:),e,v
 REAL(iwp),INTENT(OUT)::km(:,:)
 INTEGER::nod,i,j
 REAL(iwp)::aa,bb,t2,t3,t4,t5,t6,t7,t8,t10,t11,t15,t17,t18,t19,t20,t22,   &
   t23,t25,t29,t31,t33,t34,t37,t39,t40,t41,t45,t46,t48,t49,t50,t51,t52,   &
   t53,t56,t57,t58,t59,t61,t63,t66,t67,t72,t78,t82,t84,t87,t88,t92,t96,   &
   t97,t101,t102,t104,t107,t108,t110,t113,t117,t123,t126,t127,t134,t141,  &
   t150,d0=0.0_iwp,d1=1.0_iwp,d2=2.0_iwp,d4=4.0_iwp,d5=5.0_iwp,d6=6.0_iwp,&
   d7=7.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d12=12.0_iwp,d16=16.0_iwp,            &
   d17=17.0_iwp,d18=18.0_iwp,d24=24.0_iwp,d72=72.0_iwp
 nod=UBOUND(coord,1)
 SELECT CASE(nod)
 CASE(4)
   aa=coord(3,1)-coord(2,1)
   bb=coord(2,2)-coord(1,2)
   t2=-d1+d2*v
   t3=e*t2
   t4=aa**2
   t6=-d1+v
   t7=d2*e*t6
   t8=bb**2
   t10=t3*t4+t7*t8
   t11=d1/aa
   t15=d1/(d1+v)
   t17=d1/t2
   t18=d1/bb*t15*t17
   t20=t10*t11*t18/d6
   t23=e*t15*t17/d8
   t25=-e*t2*t4
   t31=(t25+e*t6*t8)*t11*t18/d6
   t37=e*(d4*v-d1)*t15*t17/d8
   t40=-t10*t11*t18/d12
   t46=(-t25-d4*e*t6*t8)*t11*t18/d12
   t48=t3*t8
   t49=t7*t4+t48
   t52=t49*t11*t18/d6
   t58=(-d4*e*t6*t4+t48)*t11*t18/d12
   t61=-t49*t11*t18/d12
   t67=(e*t6*t4-t48)*t11*t18/d6
   km(1,1)=t20
   km(1,2)=-t23
   km(1,3)=t31
   km(1,4)=t37
   km(1,5)=t40
   km(1,6)=t23
   km(1,7)=t46
   km(1,8)=-t37
   km(2,2)=t52
   km(2,3)=-t37
   km(2,4)=t58
   km(2,5)=t23
   km(2,6)=t61
   km(2,7)=t37
   km(2,8)=t67
   km(3,3)=t20
   km(3,4)=t23
   km(3,5)=t46
   km(3,6)=t37
   km(3,7)=t40
   km(3,8)=-t23
   km(4,4)=t52
   km(4,5)=-t37
   km(4,6)=t67
   km(4,7)=-t23
   km(4,8)=t61
   km(5,5)=t20
   km(5,6)=-t23
   km(5,7)=t31
   km(5,8)=t37
   km(6,6)=t52
   km(6,7)=-t37
   km(6,8)=t58
   km(7,7)=t20
   km(7,8)=t23
   km(8,8)=t52
 CASE(8)
   aa=coord(5,1)-coord(3,1)
   bb=coord(3,2)-coord(1,2)
   t2=-d1+d2*v
   t3=e*t2
   t4=aa**2
   t5=t3*t4
   t6=-d1+v
   t7=d2*e*t6
   t8=bb**2
   t11=d1/aa
   t15=d1/(d1+v)
   t17=d1/t2
   t18=d1/bb*t15*t17
   t19=d5*(t5+t7*t8)*t11*t18
   t20=t19/d18
   t22=e*t15*t17
   t23=d17/d72*t22
   t25=-d4+d8*v
   t29=-e*t6*t8
   t33=(-e*t25*t4-t29)*t11*t18/d9
   t34=d12*v
   t37=t15*t17
   t39=e*(t34-d1)*t37/d18
   t40=e*t6
   t41=t40*t8
   t45=(t5+t41)*t11*t18/d6
   t46=d4*v
   t50=e*(t46-d1)*t37/d24
   t51=d8*e*t6
   t53=-t5-t51*t8
   t56=t53*t11*t18/d18
   t57=t22/d18
   t58=t19/36
   t59=d7/d72*t22
   t61=e*(-d2+t46)
   t63=-t61*t4-t41
   t66=t63*t11*t18/d9
   t67=d4*e*t6
   t72=(t5+t67*t8)*t11*t18/d12
   t78=(t5-d16*e*t6*t8)*t11*t18/d18
   t82=e*(t34-d5)*t37/d18
   t84=t3*t8
   t87=d5*(t7*t4+t84)*t11*t18
   t88=t87/d18
   t92=-e*t2*t8
   t96=(-d16*e*t6*t4-t92)*t11*t18/d18
   t97=t67*t4
   t101=(t97+t84)*t11*t18/d12
   t102=t40*t4
   t104=-t102-t61*t8
   t107=t104*t11*t18/d9
   t108=t87/36
   t110=-t51*t4-t84
   t113=t110*t11*t18/d18
   t117=(t102+t84)*t11*t18/d6
   t123=(t102-e*t25*t8)*t11*t18/d9
   t126=-d4/d9*t63*t11*t18
   t127=d2/d9*t22
   t134=-d2/d9*t110*t11*t18
   t141=-d2/d9*t53*t11*t18
   t150=-d4/d9*t104*t11*t18
   km(1,1)=t20
   km(1,2)=-t23
   km(1,3)=t33
   km(1,4)=t39
   km(1,5)=t45
   km(1,6)=-t50
   km(1,7)=t56
   km(1,8)=t57
   km(1,9)=t58
   km(1,10)=-t59
   km(1,11)=t66
   km(1,12)=t57
   km(1,13)=t72
   km(1,14)=t50
   km(1,15)=t78
   km(1,16)=-t82
   km(2,2)=t88
   km(2,3)=-t82
   km(2,4)=t96
   km(2,5)=t50
   km(2,6)=t101
   km(2,7)=t57
   km(2,8)=t107
   km(2,9)=-t59
   km(2,10)=t108
   km(2,11)=t57
   km(2,12)=t113
   km(2,13)=-t50
   km(2,14)=t117
   km(2,15)=t39
   km(2,16)=t123
   km(3,3)=t126
   km(3,4)=d0
   km(3,5)=t33
   km(3,6)=t82
   km(3,7)=d0
   km(3,8)=t127
   km(3,9)=t66
   km(3,10)=t57
   km(3,11)=d4/d9*(t5+t29)*t11*t18
   km(3,12)=d0
   km(3,13)=t66
   km(3,14)=-t57
   km(3,15)=d0
   km(3,16)=-t127
   km(4,4)=t134
   km(4,5)=-t39
   km(4,6)=t96
   km(4,7)=t127
   km(4,8)=d0
   km(4,9)=t57
   km(4,10)=t113
   km(4,11)=d0
   km(4,12)=d2/d9*(t97+t92)*t11*t18
   km(4,13)=-t57
   km(4,14)=t113
   km(4,15)=-t127
   km(4,16)=d0
   km(5,5)=t20
   km(5,6)=t23
   km(5,7)=t78
   km(5,8)=t82
   km(5,9)=t72
   km(5,10)=-t50
   km(5,11)=t66
   km(5,12)=-t57
   km(5,13)=t58
   km(5,14)=t59
   km(5,15)=t56
   km(5,16)=-t57
   km(6,6)=t88
   km(6,7)=-t39
   km(6,8)=t123
   km(6,9)=t50
   km(6,10)=t117
   km(6,11)=-t57
   km(6,12)=t113
   km(6,13)=t59
   km(6,14)=t108
   km(6,15)=-t57
   km(6,16)=t107
   km(7,7)=t141
   km(7,8)=d0
   km(7,9)=t78
   km(7,10)=t39
   km(7,11)=d0
   km(7,12)=-t127
   km(7,13)=t56
   km(7,14)=-t57
   km(7,15)=d2/d9*(-t5+d4*e*t6*t8)*t11*t18
   km(7,16)=d0
   km(8,8)=t150
   km(8,9)=-t82
   km(8,10)=t123
   km(8,11)=-t127
   km(8,12)=d0
   km(8,13)=-t57
   km(8,14)=t107
   km(8,15)=d0
   km(8,16)=d4/d9*(-t102-t92)*t11*t18
   km(9,9)=t20
   km(9,10)=-t23
   km(9,11)=t33
   km(9,12)=t39
   km(9,13)=t45
   km(9,14)=-t50
   km(9,15)=t56
   km(9,16)=t57
   km(10,10)=t88
   km(10,11)=-t82
   km(10,12)=t96
   km(10,13)=t50
   km(10,14)=t101
   km(10,15)=t57
   km(10,16)=t107
   km(11,11)=t126
   km(11,12)=d0
   km(11,13)=t33
   km(11,14)=t82
   km(11,15)=d0
   km(11,16)=t127
   km(12,12)=t134
   km(12,13)=-t39
   km(12,14)=t96
   km(12,15)=t127
   km(12,16)=d0
   km(13,13)=t20
   km(13,14)=t23
   km(13,15)=t78
   km(13,16)=t82
   km(14,14)=t88
   km(14,15)=-t39
   km(14,16)=t123
   km(15,15)=t141
   km(15,16)=d0
   km(16,16)=t150
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
 DO i=1,nod*2
   DO j=i+1,nod*2
     km(j,i)=km(i,j)
   END DO
 END DO
RETURN
END SUBROUTINE rect_km
