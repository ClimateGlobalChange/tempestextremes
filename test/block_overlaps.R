load("~/block_r_data/stats_cold_spell.RData")

check_overlaps<-function(Alt,Alb,All,Alr,
                         Blt,Blb,Bll,Blr){
  #Note: lons must deal with periodic boundary
  #use centered lons for Atlantic
  #Check that the latitudes overlap
  if ((Alt<Blb) | (Alb>Blt)){
    ###print("lats don't overlap")
    return(FALSE)
  }
  #Check that the longitudes overlap
  else if((Alr<Bll) | (All>Blr)){
    ###print("lons don't overlap")
    return(FALSE)
  }else{
    ###print("overlap!")
    return(TRUE)
  }
}


v_count<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric())
var_tally<-data.frame(PV=0,Z=0,GHG=0)
for (t in sort(unique(df_tot$tstep))){
#for (t in 129:129){
  df_sub<-df_tot[df_tot$tstep==t,c("var","date","hr","minlat","maxlat","minlon_c","maxlon_c")]
  ###print(df_sub)
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
        ###print("swapping 1 and 2")
        vtemp<-v1
        vartemp<-varA
        v1<-v2
        varA<-varB
        v2<-vtemp
        varB<-vartemp
      }
      l1<-dim(v1)[1]
      l2<-dim(v2)[1]
      
      poss_A<-l1
      poss_B<-l1
      #print(s#printf("possible A, B: %d, %d",poss_A, poss_B))
      var_tally[1,varA]<-var_tally[1,varA]+poss_A
      var_tally[1,varB]<-var_tally[1,varB]+poss_B
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
          ###print("Min: swapping 1 and 3")
          vtemp<-v1
          vartemp<-varA
          v1<-v3
          varA<-varC
          v3<-vtemp
          varC<-vartemp
        }else{
          ###print("Min: swapping 1 and 2")
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
          ###print("Max: swapping 1 and 3")
          vtemp<-v1
          vartemp<-varA
          v1<-v3
          varA<-varC
          v3<-vtemp
          varC<-vartemp
        }else{
          ###print("Min: swapping 2 and 3")
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
      
      poss_A<-2*l1
      poss_B<-l1+l2
      poss_C<-l1+l2
      
      #print(s#printf("possible A, B, C: %d, %d, %d",poss_A, poss_B, poss_C))
      
      var_tally[1,varA]<-var_tally[1,varA]+poss_A
      var_tally[1,varB]<-var_tally[1,varB]+poss_B
      var_tally[1,varC]<-var_tally[1,varC]+poss_C
      ##print("nblobs:")
      ##print(c(l1,l2,l3))
      ##print(s##printf("%s,%s,%s",varA,varB,varC))
      #ar<-array(rep(0, l1*l2*l3), dim=c(l1,l2,l3))
    }
    ###print(s##printf("NOTE:t=%d",t))
    #now that the variables are sorted, check overlaps
    A_MAT<-matrix(0,nrow=l1,ncol=3)
    colnames(A_MAT)<-c("AB","AC","ABC")
    B_MAT<-matrix(0,nrow=l2,ncol=1)

    for (a in 1:dim(v1)[1]){
      check_AB<-FALSE
      for (b in 1:dim(v2)[1]){
        if (check_AB==FALSE){
          check_AB<-check_overlaps(v1$maxlat[a],v1$minlat[a],v1$minlon_c[a],v1$maxlon_c[a],
                                   v2$maxlat[b],v2$minlat[b],v2$minlon_c[b],v2$maxlon_c[b])
        }

        str_AB<-ifelse(check_AB==TRUE,"true","false")

        if (!is.null(v3)){
          check_BC<-FALSE
          check_AC<-FALSE

          for (c in 1:dim(v3)[1]){
            ##print(s##printf("Checking %d,%d,%d:",a,b,c))
            if (check_BC==FALSE){
              check_BC<-check_overlaps(v3$maxlat[c],v3$minlat[c],v3$minlon_c[c],v3$maxlon_c[c],
                                       v2$maxlat[b],v2$minlat[b],v2$minlon_c[b],v2$maxlon_c[b])
            }

            str_BC<-ifelse(check_BC==TRUE,"true","false")
            if (check_AC==FALSE){
              check_AC<-check_overlaps(v1$maxlat[a],v1$minlat[a],v1$minlon_c[a],v1$maxlon_c[a],
                                       v3$maxlat[c],v3$minlat[c],v3$minlon_c[c],v3$maxlon_c[c])
            }

            str_AC<-ifelse(check_AC==TRUE,"true","false")
            ##print(s##printf("AB, AC, BC: %s, %s, %s",str_AB,str_AC,str_BC))
            
            if (check_AB==TRUE & check_AC==TRUE & check_BC==TRUE){
              ##print("BREAK!")
              break
            }

            ##print(s##printf("Ending c=%d",c))
          }
          if (check_AB==TRUE){
            if (check_AC==TRUE){
              if (check_BC==TRUE){
                #ABC
                A_MAT[a,3]=1
              }
              else if (check_BC==FALSE){
                #AB and AC
                A_MAT[a,1]=1
                A_MAT[a,2]=1
              }
            }
            else if (check_AC==FALSE){
              if (check_BC==TRUE){
                #AB and BC
                A_MAT[a,1]=1
                B_MAT[b,1]=1
              }
              else if (check_BC==FALSE){
                #AB
                A_MAT[a,1]=1
              }
            }
          }else if (check_AB==FALSE){
            if (check_AC==TRUE){
              if (check_BC==TRUE){
                #AC and BC
                A_MAT[a,2]=1
                B_MAT[b,1]=1
              }
              else if (check_BC==FALSE){
                #AC
                A_MAT[a,2]=1
              }
            }
            else if (check_AC==FALSE & check_BC==TRUE){
              #BC
              B_MAT[b,1]=1
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

#Number of blocks
nPV<-dim(df_tot[df_tot$var=="PV",])[1]
nZ<-dim(df_tot[df_tot$var=="Z",])[1]
nGHG<-dim(df_tot[df_tot$var=="GHG",])[1]
nTOT<-nPV+nZ+nGHG

#Percentage of blocks
PPV<-nPV/nTOT
PZ<-nZ/nTOT
PGHG<-nGHG/nTOT

#Tally of possible blocks per overlaps is in var_tally
#Actual number of blocks per overlaps

  
nZZG<-dim(v_count[v_count$Z==1 & v_count$GHG==1 ,])[1]
nPVZ<-dim(v_count[v_count$Z==1  & v_count$PV==1,])[1]
nPVZG<-dim(v_count[v_count$GHG==1 & v_count$PV==1,])[1]
nall<-dim(v_count[v_count$Z==1 & v_count$GHG==1 & v_count$PV==1,])[1]

###print("num combo:")
###print(c(nPVZ,nPVZG,nZZG,nall))

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

###print("prob combo:")
###print(c(PZPV,PPVGHG,PZGHG,Pall))

matp<-matrix(1,nrow=3,ncol=3)

matp[1,2]<-PZPV/PZ
matp[1,3]<-PPVGHG/PGHG
matp[2,1]<-PZPV/PPV
matp[2,3]<-PZGHG/PGHG
matp[3,1]<-PPVGHG/PPV
matp[3,2]<-PZGHG/PZ



