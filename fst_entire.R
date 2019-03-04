library(ggplot2)
library(ggfortify)



#mix4entire=mix
loc=c("Oxford")
bim=read.table("/Users/Tim/Desktop/xt/1kg_entire/1kg.bim",header = F)
frqst=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/CEU.frq",header = T)
colnames(frqst)[5]="CEU"
f1=cbind(frqst[,1:2],bim[,4],frqst[,3:6])
colnames(f1)[3]="LOC"
f1=f1[-801359,]
f2=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/CHB.frq",header = T)
colnames(f2)[5]="CHB"
f3=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/YRI.frq",header = T)
colnames(f3)[5]="YRI"
f4=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/FIN.frq",header = T)
colnames(f4)[5]="FIN"
f5=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/GBR.frq",header = T)
colnames(f5)[5]="GBR"
f6=read.table("/Users/Tim/Desktop/xt/1kg_entire/standfrq/TSI.frq",header = T)
colnames(f6)[5]="TSI"
ox=read.table("/Users/Tim/Desktop/bioinfo/bb_British/frq/Oxford.frq",header = T)
colnames(ox)[5]="Oxford"



num=4
x=fst(num,loc)
frq=x[,c(2,3,6:ncol(x))]
frq=frq[order(frq[,1],frq[,2]),]

seq=5000000 #每个区间包含的kb数
index=matrix(0,22,100)
for(i in 1:22){
  loc=frq[frq[,1]==i,2]
  kb=0
  count=1
  while(kb<loc[length(loc)]){
    n=sum(loc>=kb&loc<kb+seq)
    index[i,count]=n
    kb=kb+seq
    count=count+1
  }
}

count=sum(index>10)
id=rep(0,count+1)
name=rep(" ",count+1)
tem=1
id[1]=0
for(i in 1:22){
  snp=1
  for(j in 1:100){
    if(index[i,j]>10){
      tem=tem+1
      id[tem]=index[i,j]+id[tem-1]
      #name[tem]=paste0(i,"-",snp)
      #会导致画图时第10号染色体跟在1号后面
      name[tem]=i*1000+snp
      snp=snp+1
    }
    else
      id[tem]=index[i,j]+id[tem]
  }
}
name=name[-1]
size=length(name)
st=frq[,3:(num+2)]
co=frq[,(num+3)]
z=data.frame(Chr=1:(num*size),Origin=1:(num*size),Distance=1:(num*size))
z[,1]=name
for(i in 1:num){
  for(j in 1:size){
    if(j==1){
      si=co[1:id[2]]
      s=st[1:id[2],]
    }
    else{
      si=co[(id[j]+1):id[j+1]]
      s=st[(id[j]+1):id[j+1],]
    }
    #w1=gn[r[i]]/(gn[r[i]]+gn[j])
    #w2=1-w1
    w1=w2=0.5
    pbar=w1*s[,i]+w2*si
    sum=w1*(s[,i]-pbar)^2+w2*(si-pbar)^2
    pbarnon=pbar[pbar!=0]
    sumnon=sum[pbar!=0]
    sumnon=sumnon/(2*pbarnon*(1-pbarnon))
    z[(i-1)*size+j,2]=colnames(s)[i]
    z[(i-1)*size+j,3]=mean(sumnon)
  }
  
}

z$Chr=factor(z$Chr,levels = name)
#柱状图默认按照阿斯克码排列，所以先要指定按染色体顺序排
r=z
r[,3]=1/r[,3]


lab=paste0("Chr",1:22)
brk=seq(1001,22001,1000)
ggplot(r,aes(x=Chr,y=Distance,fill=Origin)) + geom_bar(stat='identity',position="fill") +
  scale_x_discrete(breaks=brk,labels=lab) + theme(axis.text.x=element_text(size=8))



fmerge=function(m,f,AT=FALSE){
  f=f[,c(-1,-6)]
  mix=merge(m,f,by.x = 'SNP', by.y = 'SNP', all = FALSE)
  #judge 5 simulations
  a=mix[,c(4,5,ncol(mix)-2,ncol(mix)-1)]
  #without number of individual
  a=apply(a,2,ordered,levels = c('A', 'C', 'G', 'T'),labels = c('1', '2', '3','4'))
  a=apply(a,2,as.numeric)
  sim=(a[,1]==a[,3])+(a[,1]==a[,4])+(a[,2]==a[,3])+(a[,2]==a[,4])
  sim[sim!=0&sim!=2]=9
  sim[sim==2&a[,1]==a[,3]]=1
  sim[sim==0&(a[,1]+a[,3])==5]=3
  sim[sim==0&(a[,1]+a[,4])==5]=4
  sim[sim==0]=9
  mix[sim==2|sim==4,ncol(mix)]=1-mix[sim==2|sim==4,ncol(mix)]
  #remove AT/TA and CG/GC
  if(AT)  
  {
    sim[(a[,1]+a[,2])==5]=9
  }
  mix=mix[sim!=9,]
  mix=mix[,-c(ncol(mix)-2,ncol(mix)-1)]
  return(mix)
}

fst=function(n,loc,AT=FALSE){
  mix1=f1[,-7]
  # if(length(loc)==1){
  #   cof=paste0("/Users/Tim/Desktop/bioinfo/bb_British/frq/",loc,".frq")
  #   co=read.table(cof,header = T)
  #   colnames(co)[5]=as.character(loc)
  # }
  co=ox
  if(n==3){
    mix2=fmerge(mix1,f2)
    mix3=fmerge(mix2,f3)
    mixco=fmerge(mix3,co)
  }
  if(n==4){
    mix4=fmerge(mix1,f4)
    mix5=fmerge(mix4,f5)
    mix6=fmerge(mix5,f6)
    mixco=fmerge(mix6,co)
  }
  return(mixco)
}



####画每个cohort平均fst距离
loc=read.table("/Users/Tim/Desktop/bioinfo/ukb/location.txt",header = F)
x=flip(4)
frq=x[,c(2,3,6:ncol(x))]
frq=frq[order(frq[,1],frq[,2]),]
z=fstcal(frq,4)
frq[is.na(frq)]=0.3  #missing data
r=z
r[,3]=1/r[,3]
ggplot(r,aes(x=Cohort,y=Distance,fill=Origin)) + geom_bar(stat='identity',position="fill") + 
  theme(axis.text.x=element_text(size=12,hjust=1,vjust=0.5,angle=90))

# ggplot(r,aes(x=Chr,y=Distance,fill=Stand)) + geom_bar(stat='identity',position="fill") + 
#   theme(axis.text.x=element_text(size=3,hjust=1,vjust=0.5,angle=90))
# 
# ggplot(r,aes(x=Chr,y=Distance,fill=Stand)) + geom_bar(stat='identity',position="fill") +
#   theme(axis.text.x=element_blank())


flip=function(n,AT=FALSE){
  if(n==3) {
    c=c("/Users/Tim/Desktop/xt/1kg_entire/standfrq/CHB.frq","/Users/Tim/Desktop/xt/1kg_entire/standfrq/YRI.frq")
    l=c("CHB","YRI",loc)
  }
  if(n==4) {
    c=c("/Users/Tim/Desktop/xt/1kg_entire/standfrq/FIN.frq","/Users/Tim/Desktop/xt/1kg_entire/standfrq/GBR.frq","/Users/Tim/Desktop/xt/1kg_entire/standfrq/TSI.frq")
    l=c("FIN","GBR","TSI",as.character(loc[,1]))
  }
  cof=paste0("/Users/Tim/Desktop/bioinfo/bb_British/frq/",as.character(loc[,1]),".frq")
  list=c(c,cof)
  mix=f1[,-7]
  #cohort人数暂时没用
  #colnames(mix)[6]="CEU.N"
  for(i in 1:length(list)){
    add=read.table(list[i],header = T)
    colnames(add)[5]=as.character(l[i])
    #colnames(add)[6]=paste0(l[i],".N")
    mix=fmerge(mix,add)
  }
  return(mix)
}

fstcal=function(frq,num){
  s=frq[,3:(num+2)]
  n=ncol(s)
  si=frq[,(num+3):ncol(frq)]
  size=ncol(si)
  z=data.frame(Cohort=1:(n*size),Origin=1:(n*size),Distance=1:(n*size))
  z[,1]=colnames(si)
  for(i in 1:n){
    for(j in 1:size){
      #w1=gn[r[i]]/(gn[r[i]]+gn[j])
      #w2=1-w1
      w1=w2=0.5
      pbar=w1*s[,i]+w2*si[,j]
      sum=w1*(s[,i]-pbar)^2+w2*(si[,j]-pbar)^2
      pbarnon=pbar[pbar!=0]
      sumnon=sum[pbar!=0]
      sumnon=sumnon/(2*pbarnon*(1-pbarnon))
      z[(i-1)*size+j,2]=colnames(s)[i]
      z[(i-1)*size+j,3]=mean(sumnon)
    }
  }
  return(z)
}

