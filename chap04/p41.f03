PROGRAM p41
!-------------------------------------------------------------------------
! Program 4.1 One dimensional analysis of axially loaded elastic rods
!             using 2-node rod elements.
! 注释者：wangyb19     时间：2022/11/02
!-------------------------------------------------------------------------
 USE main
 USE geom
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndof=2,nels,neq,nlen,nod=2, &
   nodof=1,nn,nprops=1,np_types,nr
!-------------------------------------------------------------------------
! fixed_freedoms: 约束位移的自由度数（位移边界条件，区别于字面意思fixed）
! i,iel: 计数器
! k: 固定位移的节点编号
! loaded_nodes: 力边界条件的节点数目
! ndof=2:每个单元的自由度数 
! nels: 整个模型的单元数目
! neq: 模型自由度数
! nlen: 文件名长度 
! nod=2: 每个单元的节点数
! nodof=1: 每个单元有一个自由度
! nn: 总节点数
! nprops=1: 本构参数数目，本模型只有一个：杨氏模量
! np_types: 本构类型数目
!-------------------------------------------------------------------------
 REAL(iwp)::penalty=1.0e20_iwp,zero=0.0_iwp
 CHARACTER(LEN=15)::argv
!-------------------------------------------------------------------------
! argv: 输入文件名
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),kdiag(:),nf(:,:),no(:),      &
   node(:),num(:) 
 REAL(iwp),ALLOCATABLE::action(:),eld(:),ell(:),km(:,:),kv(:),loads(:),   &
   prop(:,:),value(:)
!-----------------------input and initialisation--------------------------
! CALL getname(argv,nlen)
 argv="p41_1"
 nlen=5
 OPEN(10,FILE=argv(1:nlen)//'.dat')
 OPEN(11,FILE=argv(1:nlen)//'.res')
 READ(10,*)nels,np_types
 nn=nels+1 ! 根据输入的单元数目计算得到节点数目
 ALLOCATE(g(ndof),num(nod),nf(nodof,nn),etype(nels),ell(nels),eld(ndof),  &
   km(ndof,ndof),action(ndof),g_g(ndof,nels),prop(nprops,np_types))
!-------------------------------------------------------------------------
! g: 一维数组 存储单元局部的节点编号情况
! num: 不考虑固定约束时每个单元的节点编号信息
! nf: 模型整体的节点编号信息（固定约束赋零值）
! etype: 存储每个单元的性质编号
! ell: 存储每个单元的长度信息
! eld: 存储每个单元的位移信息
! km: 单元刚度矩阵
! action: 存储每个单元的应力信息
! g_g: 存储所有单元的节点编号情况
! prop: 存储所有单元的物理模量
!-------------------------------------------------------------------------
 READ(10,*)prop
 etype=1
 IF(np_types>1)READ(10,*)etype
 READ(10,*)ell
 nf=1
 READ(10,*)nr,(k,nf(:,k),i=1,nr)
 CALL formnf(nf) ! 给单元编号，固定约束赋零值
 neq=MAXVAL(nf)  ! 根据上一步的编号情况得到系统的总自由度数
 ALLOCATE(kdiag(neq),loads(0:neq))
 kdiag=0
!-----------------------loop the elements to find global arrays sizes-----
 elements_1: DO iel=1,nels
   num=(/iel,iel+1/) ! 单元默认的节点编号（未考虑固定约束）
   CALL num_to_g(num,nf,g) ! 考虑固定约束后生成每个单元的节点编号信息
   g_g(:,iel)=g
   CALL fkdiag(kdiag,g) ! skyline 刚度矩阵每一行的非零元素个数（考虑对称性）
 END DO elements_1
 DO i=2,neq
   kdiag(i)=kdiag(i)+kdiag(i-1)
 END DO
 ALLOCATE(kv(kdiag(neq))) !以列向量的形式储存整体刚度矩阵的所有非零元素
 WRITE(11,'(2(A,I5))')                                                    &
   " There are",neq," equations and the skyline storage is",kdiag(neq)
!-----------------------global stiffness matrix assembly------------------
 kv=zero
 elements_2: DO iel=1,nels
   CALL rod_km(km,prop(1,etype(iel)),ell(iel)) ! 生成单元刚度矩阵
   g=g_g(:,iel)   ! 从全局节点编号情况里摘取局部节点编号情况
   CALL fsparv(kv,km,g,kdiag) ! 以向量形式组装整体刚度矩阵
 END DO elements_2
!-----------------------read loads and/or displacements-------------------
 loads=zero
 READ(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)!读受载节点情况
 READ(10,*)fixed_freedoms
 IF(fixed_freedoms/=0)THEN
   ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),value(fixed_freedoms))
   !----------------------------------------------------------------------
   ! node
   ! no: 存储施加位移边界条件的节点的编号信息
   !----------------------------------------------------------------------
   READ(10,*)(node(i),value(i),i=1,fixed_freedoms)
   DO i=1,fixed_freedoms
     no(i)=nf(1,node(i))
   END DO
   kv(kdiag(no))=kv(kdiag(no))+penalty ! 对刚度元素处理以引入位移边界条件
   loads(no)=kv(kdiag(no))*value 
 END IF
!-----------------------equation solution -------------------------------- 
 CALL sparin(kv,kdiag)
 CALL spabac(kv,loads,kdiag)
 loads(0)=zero
 WRITE(11,'(/A)')"  Node   Disp"
 DO k=1,nn
   WRITE(11,'(I5,2E12.4)')k,loads(nf(:,k))
 END DO
!-----------------------retrieve element end actions----------------------
 WRITE(11,'(/A)')" Element Actions"
 elements_3: DO iel=1,nels
   CALL rod_km(km,prop(1,etype(iel)),ell(iel))
   g=g_g(:,iel)
   eld=loads(g)
   action=MATMUL(km,eld)
   WRITE(11,'(I5,2E12.4)')iel,action
 END DO elements_3
STOP
END PROGRAM p41
