#load("~/block_r_data/stats_cold_spell.RData")

check_overlaps<-function(Alt,Alb,All,Alr,
                         Blt,Blb,Bll,Blr){
  #Note: lons must deal with periodic boundary
  #use centered lons for Atlantic
  #Check that the latitudes overlap
  if ((Alt<Blb) | (Alb>Blt)){
    ##print("lats don't overlap")
    return(FALSE)
  }
  #Check that the longitudes overlap
  else if((Alr<Bll) | (All>Blr)){
    ##print("lons don't overlap")
    return(FALSE)
  }else{
    ##print("overlap!")
    return(TRUE)
  }
}


v_count<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric())

for (t in sort(unique(df_tot$tstep))){
  #for (t in 271:271){
  df_sub<-df_tot[df_tot$tstep==t,c("var","date","hr","minlat","maxlat","minlon_c","maxlon_c")]
  ##print(df_sub)
  namevec<-unique(df_sub$var)
  if (length(namevec)>1){
    v1<-df_sub[df_sub$var==namevec[1],]
    l1<-dim(v1)[1]
    varA<-as.character(namevec[1])
    v2<-df_sub[df_sub$var==namevec[2],]
    l2<-dim(v2)[1]
    varB<-as.character(namevec[2])
    
    if (length(namevec)==2){
      v3<-NULL
      varC<-NULL
      #Compare the lengths of v1 and v2
      minlen<-min(l1,l2)
      #if v2<v1, swap, else keep the same
      if (minlen!=l1){
        ##print("swapping 1 and 2")
        vtemp<-v1
        vartemp<-varA
        v1<-v2
        varA<-varB
        v2<-vtemp
        varB<-vartemp
      }
      l1<-dim(v1)[1]
      l2<-dim(v2)[1]
      #ar<-array(rep(0, l1*l2), dim=c(l1,l2,1))
    }
    else{
      v3<-df_sub[df_sub$var==namevec[3],]
      l3<-dim(v3)[1]
      varC<-as.character(namevec[3])
      
      minlen<-min(l1,l2,l3)
      if (minlen!=l1){
        #is l2 or l3 smaller?
        m2<-min(l2,l3)
        #l3 is smallest, swap 1 and 3
        if (m2!=l2){
          ##print("Min: swapping 1 and 3")
          vtemp<-v1
          vartemp<-varA
          v1<-v3
          varA<-varC
          v3<-vtemp
          varC<-vartemp
        }else{
          ##print("Min: swapping 1 and 2")
          vtemp<-v1
          vartemp<-varA
          v1<-v2
          varA<-varB
          v2<-vtemp
          varB<-vartemp
        }
      }
      l1<-dim(v1)[1]
      l2<-dim(v2)[1]
      l3<-dim(v3)[1]
      maxlen<-max(l1,l2,l3)
      if (maxlen!=l3){
        #is l1 or l2 bigger?
        m2<-max(l1,l2)
        #l1 is biggest, swap 1 and 3
        if (m2!=l2){
          ##print("Max: swapping 1 and 3")
          vtemp<-v1
          vartemp<-varA
          v1<-v3
          varA<-varC
          v3<-vtemp
          varC<-vartemp
        }else{
          ##print("Min: swapping 2 and 3")
          vtemp<-v3
          vartemp<-varC
          v3<-v2
          varC<-varB
          v2<-vtemp
          varB<-vartemp
        }
      }
      l1<-dim(v1)[1]
      l2<-dim(v2)[1]
      l3<-dim(v3)[1]
      #print("nblobs:")
      #print(c(l1,l2,l3))
      #ar<-array(rep(0, l1*l2*l3), dim=c(l1,l2,l3))
    }
    #print(s#printf("NOTE:t=%d",t))
    #now that the variables are sorted, check overlaps
    A_MAT<-matrix(0,nrow=l1,ncol=3)
    colnames(A_MAT)<-c("AB","AC","ABC")
    B_MAT<-matrix(0,nrow=l2,ncol=1)
    for (a in 1:dim(v1)[1]){
      for (b in 1:dim(v2)[1]){
        check_AB<-FALSE
        check_AB<-check_overlaps(v1$maxlat[a],v1$minlat[a],v1$minlon_c[a],v1$maxlon_c[a],
                                 v2$maxlat[b],v2$minlat[b],v2$minlon_c[b],v2$maxlon_c[b])
        
        if (!is.null(v3)){
          for (c in 1:dim(v3)[1]){
            check_BC<-FALSE
            check_BC<-check_overlaps(v3$maxlat[c],v3$minlat[c],v3$minlon_c[c],v3$maxlon_c[c],
                                     v2$maxlat[b],v2$minlat[b],v2$minlon_c[b],v2$maxlon_c[b])
            check_AC<-FALSE
            check_AC<-check_overlaps(v1$maxlat[a],v1$minlat[a],v1$minlon_c[a],v1$maxlon_c[a],
                                     v3$maxlat[c],v3$minlat[c],v3$minlon_c[c],v3$maxlon_c[c])
            if (check_AB==TRUE){
              if (check_BC==TRUE){
                #ABC
                if (check_AC==TRUE){
                  A_MAT[a,3]=1
                }
                else{
                  #AB and BC
                  A_MAT[a,1]=1
                  B_MAT[b,1]=1
                }
              }else{
                if (check_AC==TRUE){
                  #AB and AC
                  A_MAT[a,1]=1
                  A_MAT[a,2]=1
                }else{
                  #AB
                  A_MAT[a,1]=1
                }
              }
            }else{
              if (check_BC==TRUE){
                if (check_AC==TRUE){
                  #AC and BC
                  A_MAT[a,2]=1
                  B_MAT[b,1]=1
                }
                else{
                  #BC
                  B_MAT[b,1]=1
                }
              }
              else{
                #AC
                if (check_AC==TRUE){
                  A_MAT[a,2]=1
                }
              }
            }
          }
        }
        #Two variables
        else if (check_AB==TRUE){
          A_MAT[a,1]=1
        }
      }
    }
    #Now to add a row for each of the counts
    #AB
    for (ab in 1:sum(A_MAT[,1])){
      v_add<-data.frame(PV=0,Z=0,GHG=0)
      v_add[1,varA]=1
      v_add[1,varB]=1
      v_count=rbind(v_count,v_add)
    }
    if (!is.null(v3)){
      #AC
      if (sum(A_MAT[,2])>0){
        for (ac in 1:sum(A_MAT[,2])){
          v_add<-data.frame(PV=0,Z=0,GHG=0)
          v_add[1,varA]=1
          v_add[1,varC]=1
          v_count=rbind(v_count,v_add)
        }
      }
      #ABC
      if (sum(A_MAT[,3])>0){
        for (abc in 1:sum(A_MAT[,3])){
          v_add<-data.frame(PV=1,Z=1,GHG=1)
          v_count=rbind(v_count,v_add)
        }
      }
      #BC
      if (sum(B_MAT[,1])>0){
        for (bc in 1:sum(B_MAT[,1])){
          v_add<-data.frame(PV=0,Z=0,GHG=0)
          v_add[1,varB]=1
          v_add[1,varC]=1
          v_count=rbind(v_count,v_add)
        }
      }
    }
  }
  
}

#Now process the array
nPV<-dim(df_tot[df_tot$var=="PV",])[1]
nZ<-dim(df_tot[df_tot$var=="Z",])[1]
nGHG<-dim(df_tot[df_tot$var=="GHG",])[1]
nTOT<-nPV+nZ+nGHG

##print("num each:")
##print(c(nPV,nZ,nGHG,nTOT))

PPV<-nPV/nTOT
PZ<-nZ/nTOT
PGHG<-nGHG/nTOT


# PPV<-nPV/nTimeSteps
# PZ<-nZ/nTimeSteps
# PGHG<-nGHG/nTimeSteps

##print("prob each:")
##print(c(PPV,PZ,PGHG))

nZZG<-dim(v_count[v_count$Z==1 & v_count$GHG==1 & v_count$PV==0,])[1]
nPVZ<-dim(v_count[v_count$Z==1 & v_count$GHG==0 & v_count$PV==1,])[1]
nPVZG<-dim(v_count[v_count$Z==0 & v_count$GHG==1 & v_count$PV==1,])[1]
nall<-dim(v_count[v_count$Z==1 & v_count$GHG==1 & v_count$PV==1,])[1]

##print("num combo:")
##print(c(nPVZ,nPVZG,nZZG,nall))

#PV and Z
PZPV<-nPVZ/nTOT
#PV and GHG
PPVGHG<-nPVZG/nTOT
#Z and GHG
PZGHG<-nZZG/nTOT
#all
Pall<-nall/nTOT



# #PV and Z
# PZPV<-nPVZ/nTimeSteps
# #PV and GHG
# PPVGHG<-nPVZG/nTimeSteps
# #Z and GHG
# PZGHG<-nZZG/nTimeSteps
# #all
# Pall<-nall/nTimeSteps

##print("prob combo:")
##print(c(PZPV,PPVGHG,PZGHG,Pall))

matp<-matrix(1,nrow=3,ncol=3)

matp[1,2]<-PZPV/PZ
matp[1,3]<-PPVGHG/PGHG
matp[2,1]<-PZPV/PPV
matp[2,3]<-PZGHG/PGHG
matp[3,1]<-PPVGHG/PPV
matp[3,2]<-PZGHG/PZ



