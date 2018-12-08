#load("~/block_r_data/stats_cold_spell.RData")
library(knitr)
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

#df_tot<-df_tot_nostitch
df_tot$date_hr<-sprintf("%s_%02d",df_tot$date,df_tot$hr)
v_type<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric())
for (t in sort(unique(df_tot$date_hr))){

  df_sub<-df_tot[df_tot$date_hr==t,c("var","date","hr",
                                     "minlat","maxlat",
                                     "minlon","maxlon",
                                     "minlon_c","maxlon_c")]
  if (sector=="NA" | sector == "SA"){
    minlon_check<- "minlon_c"
    maxlon_check<-"maxlon_c"
  }else{
    minlon_check<- "minlon"
    maxlon_check<-"maxlon"
  }
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

    }
    else{
      v3<-df_sub[df_sub$var==namevec[3],]
      l3<-dim(v3)[1]
      varC<-as.character(namevec[3])
      
    }
    #holds the pairs
    AB<-data.frame(A=numeric(),B=numeric())
    AC<-NULL
    BC<-NULL
    ABC<-NULL
    #New frame to hold copy of A,B,C
    AllBlobs<-data.frame(type=character(),num=numeric(),stringsAsFactors = FALSE)
    AllBlobs<-rbind(
      data.frame(type=rep("A",l1),num=seq(1,l1)),
      data.frame(type=rep("B",l2),num=seq(1,l2))
    )
    X<-data.frame(type=character(),num=numeric(), stringsAsFactors = FALSE)

    if (!is.null(v3)){
      AC<-data.frame(A=numeric(),C=numeric())
      BC<-data.frame(B=numeric(),C=numeric())
      ABC<-data.frame(A=numeric(),B=numeric(),C=numeric())
      
      AllBlobs<-rbind(
        AllBlobs,
        data.frame(type=rep("C",l3),num=seq(1,l3))
      )
    }

    #check AB
    for (a in 1:l1){
      for (b in 1:l2){
        check_AB<-check_overlaps(v1[a,"maxlat"],v1[a,"minlat"],v1[a,minlon_check],v1[a,maxlon_check],
                                 v2[b,"maxlat"],v2[b,"minlat"],v2[b,minlon_check],v2[b,maxlon_check])
        if(check_AB==TRUE){
          #Add pair to AB
          AB=rbind(AB,data.frame(A=a,B=b))
        }
      }
    }
    #check AC and BC
    if (!is.null(v3)){
      #Check AC
      for (a in 1:l1){
        for (c in 1:l3){
          check_AC<-check_overlaps(v1[a,"maxlat"],v1[a,"minlat"],v1[a,minlon_check],v1[a,maxlon_check],
                                   v3[c,"maxlat"],v3[c,"minlat"],v3[c,minlon_check],v3[c,maxlon_check])
          if (check_AC==TRUE){
            AC=rbind(AC,data.frame(A=a,C=c))
          }
        }
      }
      #check BC
      for (b in 1:l2){
        for (c in 1:l3){
          check_BC<-check_overlaps(v3[c,"maxlat"],v3[c,"minlat"],v3[c,minlon_check],v3[c,maxlon_check],
                                   v2[b,"maxlat"],v2[b,"minlat"],v2[b,minlon_check],v2[b,maxlon_check])
          if (check_BC==TRUE){
            BC=rbind(BC,data.frame(B=b,C=c))
          }
        }
      }
      
    }
    #print("BEFORE X:")
    #print(AB)
    #print(AC)
    #print(BC)
    #print("Entering while loop")
    while(dim(AllBlobs)[1]>0){
      #Generate random integer to pull out and check
      randint<-as.integer(runif(1,1,dim(AllBlobs)[1]))
      blob_check1<-AllBlobs[randint,]

      X<-rbind(X,blob_check1)
      #Pop out the blob that will be checked
      AllBlobs<-AllBlobs[-randint,]
      blob1<-as.character(blob_check1$type[1])
      b1num<-blob_check1$num[1]
      #Check against other two types
      if (blob1=="A"){
        check2<-AB
        blob2<-"B"
        check3<-BC
        blob3<-"C"
        check4<-AC
      }
      else if (blob1=="B"){
        check2<-AB
        blob2<-"A"
        check3<-AC
        blob3<-"C"
        check4<-BC
      }
      else{
        check2<-AC
        blob2<-"A"
        check3<-AB
        blob3<-"B"
        check4<-BC
      }
      #We're doing 3 checks
      #First check is to find the row where the blob matches in check2
      check2_rows<-which(check2[,blob1]==b1num)
      if (length(check2_rows)>0){
        check2_vals<-check2[check2_rows,blob2]
        #Second check is where blob matches in check3
        for (z in 1:length(check2_vals)){
          #Blob of 2nd type
          row_check2<-which(AllBlobs$type==blob2 & AllBlobs$num==check2_vals[z])
          blob_check2<-AllBlobs[row_check2,]
          X<-rbind(X,blob_check2)
          AllBlobs<-AllBlobs[-row_check2,]
          b2num<-blob_check2$num[1]
          #Check against the third type
          check3_rows<-which(check3[,blob2]==b2num)
          if (length(check3_rows)>0){
            check3_vals<-check3[check3_rows,blob3]
            for (y in 1:length(check3_vals)){
              row_check3<-which(AllBlobs$type==blob3 & AllBlobs$num==check3_vals[z])
              blob3_pop<-AllBlobs[row_check3,]
              X<-rbind(X,blob3_pop)
              AllBlobs<-AllBlobs[-row_check3,]
            }
          }
        }
      }
      #Also need to run check4 to do comparison of 1 and 3
      check4_rows<-which(check4[,blob1]==b1num)
      if (length(check4_rows)>0){
        check4_vals<-check4[check4_rows,blob3]
        for (w in 1:length(check4_vals)){
          row_check4<-which(AllBlobs$type==blob3 & AllBlobs$num==check4_vals[w])
          if (length(row_check4)>0){
            blobs4_pop<-AllBlobs[row_check4,]
            X<-rbind(X,blobs4_pop)
            AllBlobs<-AllBlobs[-row_check4,]
          }
        }
      }
      type_blobs<-unique(X$type)
      if (length(type_blobs)>1){
        #print("End result:")
        #print(X)
        v_add<-data.frame(PV=0,Z=0,GHG=0)
        for (var in type_blobs){
          v_char<-as.character(var)
          varname<-ifelse(v_char=="A",varA,ifelse(v_char=="B",varB,varC))
          v_add[1,varname]=1
        }

        v_type<-rbind(v_type,v_add)
      }
      #print("Current nrows of AllBlobs:")
      #print(dim(AllBlobs)[1])
      X<-data.frame(type=character(),num=numeric(), stringsAsFactors = FALSE)

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

  
nZZG<-dim(v_type[v_type$Z==1 & v_type$GHG==1 ,])[1]
nPVZ<-dim(v_type[v_type$Z==1  & v_type$PV==1,])[1]
nPVZG<-dim(v_type[v_type$GHG==1 & v_type$PV==1,])[1]
nall<-dim(v_type[v_type$Z==1 & v_type$GHG==1 & v_type$PV==1,])[1]
tot_overlaps<-dim(v_type)[1]

nPVoverlap<-sum(v_type$PV)
nZoverlap<-sum(v_type$Z)
nGHGoverlap<-sum(v_type$GHG)

varvec<-c("PV*","Z*","ZG")
varvec_f<-c("PV&ast;","Z&ast;","ZG")

matnum<-matrix(1,nrow=2,ncol=4)
colnames(matnum)<-c(varvec,"tot")
rownames(matnum)<-c("number","percent")
matnum[1,1]<-nPV
matnum[1,2]<-nZ
matnum[1,3]<-nGHG
matnum[1,4]<-nTOT
matnum[2,]<-round(matnum[1,]/nTOT,4)

matnum[1,]<-format(matnum[1,],nsmall=0)

cat("#### Number of instantaneous blocks by type\n")
print(kable(matnum,format='markdown',align='c'))

matint<-matrix(1,nrow=2,ncol=5)
colnames(matint)<-c("PV&ast;$\\cap$Z&ast;","PV&ast;$\\cap$ZG","Z&ast;$\\cap$ZG","all","total")
rownames(matint)<-c("number","percent")
matint[1,1]<-nPVZ
matint[1,2]<-nPVZG
matint[1,3]<-nZZG
matint[1,4]<-nall
matint[1,5]<-tot_overlaps
matint[2,]<-round(matint[1,]/tot_overlaps,4)
matint[1,]<-format(matint[1,],nsmall=0)
cat("#### Number of overlapping blocks\n")
print(kable(matint,format='markdown',align='c'))


cat("**Percent of blocks that have some sort of overlap:**",tot_overlaps/nTOT,"\n")

matp<-matrix(1,nrow=3,ncol=3)

matp[1,2]<-nPVZ/nZ
matp[1,3]<-nPVZG/nGHG
matp[2,1]<-nPVZ/nPV
matp[2,3]<-nZZG/nGHG
matp[3,1]<-nPVZG/nPV
matp[3,2]<-nZZG/nZ

mats<-matrix(1,nrow=3,ncol=3)
rownames(mats)<-varvec
colnames(mats)<-varvec


for (a in 1:3){
  v1<-varvec_f[a]
  for (b in 1:3){
    v2<-varvec_f[b]
    mats[a,b]<-sprintf("P(%s|%s)=%.4f",v1,v2,matp[a,b])
  }
}
cat("\n")
cat("#### Probability of block type given presence of other block type\n")
print(kable(mats,format='markdown',align='c'))
