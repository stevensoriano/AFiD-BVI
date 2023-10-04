      subroutine allocate_trigeo
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      implicit none 
      integer :: merr
      !-------------------------------------------------

      allocate(n_edge_of_vert(maxnv,Nparticle))
      allocate(vert_of_edge(2,maxne,Nparticle))
      allocate(face_of_edge(2,maxne,Nparticle))
      allocate(vert_of_face(3,maxnf,Nparticle))
      allocate(edge_of_face(3,maxnf,Nparticle))
      allocate(vert_of_vert(max_n_edge_of_vert,maxnv,Nparticle))
      allocate(edge_of_vert(max_n_edge_of_vert,maxnv,Nparticle))
      allocate(v1234(4,maxne,Nparticle))

      allocate(tri_ver(9,maxnf,Nparticle))
      allocate(vel_tri(3,maxnf,Nparticle),acc_tri(3,maxnf,Nparticle))
      allocate(tri_bar(3,maxnf,Nparticle),tri_nor(3,maxnf,Nparticle))
      allocate(bboxind(3,2,Nparticle))
      allocate(dum_for(maxnf,Nparticle))
      allocate(pind(6,maxnf,Nparticle))
      allocate(mvol(maxnf,Nparticle),press_face(maxnf,Nparticle))
      allocate(vforc_face(maxnf,Nparticle))

      allocate(theta0(maxne,Nparticle),theta(maxne,Nparticle))
      allocate(dist0(maxne,Nparticle),dist(maxne,Nparticle))
      allocate(sur0(maxnf,Nparticle),sur(maxnf,Nparticle))
      allocate(dismax(maxnf,Nparticle))

      allocate(Volume0(Nparticle),Volume(Nparticle))
      allocate(Surface0(Nparticle),Surface(Nparticle))
      allocate(shwtx(7),shwty(7),shwtz1(7),shwtz2(7))

      allocate(xyz0(3,maxnv,Nparticle),xyz(3,maxnv,Nparticle))
      allocate(xyzv(3,maxnv,Nparticle),xyza(3,maxnv,Nparticle))
      allocate(xyzv1(3,maxnv,Nparticle),xyza0(3,maxnv,Nparticle))
      
      allocate(fpxyz(3,maxnv,Nparticle),fexyz(3,maxnv,Nparticle))
      allocate(fbxyz(3,maxnv,Nparticle),fvxyz(3,maxnv,Nparticle))
      allocate(fatxyz(3,maxnv,Nparticle),falxyz(3,maxnv,Nparticle))
      allocate(fsxyz(3,maxnv,Nparticle),fxyz(3,maxnv,Nparticle))

      allocate(sca_nod1(maxnf,Nparticle))
      allocate(sca_nod2(maxnf,Nparticle))
      allocate(sca_nod3(maxnf,Nparticle))

      return
      end 
