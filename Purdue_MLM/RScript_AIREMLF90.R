
# AIREMLF90 available at http://nce.ads.uga.edu/html/projects/programs/

# Function to create inputs for remlf90
# gen: numeric matrix with row names as ID - coded as 012
# data: data.frame with a column called ID
# form: formula

ram = function(form, data, gen, Indep_Res=FALSE){
  
  # Part 1
  ter = terms(form)
  labs = labels(ter)
  labs = labs[labs!='ID']
  f = length(labs)
  classes = levs = rep(1,f)
  for(i in 1:f){
    classes[i]=class(data[[labs[i]]])
    if(classes[i]=='factor') levs[i] = length(levels(data[[labs[i]]]))
  }
  classes = ifelse(classes=='factor','cross','cov')
  # Part 2
  lines0 = levels(data$ID)
  lines1 = rownames(gen)
  lines2 = 1:length(lines0)
  GinO = lines1%in%lines0
  OinG = lines0%in%lines1
  lines3 = lines2[OinG] # goes to geno
  cat('Building genotypic and pedigree files\n')
  gen[is.na(gen)]=5
  c0125 = apply(gen[GinO,],1,paste,collapse='')
  spaces = sprintf("%-8s", lines3)
  geno = data.frame(spaces,c0125)
  write.table(geno,file='geno',sep=' ',row.names=F,col.names=F,quote=F)
  xref = data.frame(lines2,lines2)
  write.table(xref,file='geno_XrefID',sep=' ',row.names=F,col.names=F,quote=F)
  ped = data.frame(lines2,rep(0,length(lines0)),rep(0,length(lines0)))
  write.table(ped,file='ped',sep=' ',row.names=F,col.names=F,quote=F)
  index = data.frame(lines0,lines2)
  write.table(index,file='name_index',sep=' ',row.names=F,col.names=F,quote=F)
  # Part 3
  cat('Building data and parameter files\n')
  data2 = data.frame(rep(1,nrow(data)))
  for(i in 1:f) data2 = cbind(data2,as.numeric(data[[labs[i]]]))
  Ys = deparse(ter[[2]]);Ys=unlist(strsplit(Ys,' '));Ys=Ys[Ys!='|']
  tr = length(Ys)
  data2 = data.frame(data2,as.numeric(data$ID))
  for(i in 1:tr) data2 = data.frame(data2,as.numeric(data[[Ys[i]]]))
  data2[is.na(data2)] = 0
  write.table(data2,file='data',sep=' ',row.names=F,col.names=F,quote=F)
  # Part 4
  classes = c('cov',classes,'cross'); levs = c(1,levs,length(levels(data$ID)))
  para = c()
  para = c(para,
           'DATAFILE
data
NUMBER_OF_TRAITS',
           tr,
           'NUMBER_OF_EFFECTS',
           (f+2),
           'OBSERVATION(S)',
           paste((f+3):(f+2+tr),collapse = ' '),
           'WEIGHT(S)
           
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT [EFFECT NESTED]')
  for(i in 1:(f+2)) para=c(para,paste(c(paste(rep(i,tr),collapse=' '),levs[i],classes[i]),collapse=' '))
  para=c(para,'RANDOM_RESIDUAL VALUES')
  for(i in 1:tr){
    if(Indep_Res){
      Row = rep(0,tr); Row[i]=1; Row=paste(Row,collapse=' ')
    }else{
      Row = rep(0.1,tr); Row[i]=1.0; Row=paste(Row,collapse=' ')
    }
    para = c(para,Row)}
  for(i in which(classes=='cross')){
    u = c('RANDOM_GROUP',i,
          'RANDOM_TYPE',
          ifelse(i!=(f+2),'diagonal','add_animal'),
          'FILE',
          ifelse(i!=(f+2),'','ped'),
          '(CO)VARIANCES')
    for(j in 1:tr){
      Row = rep( ifelse(i!=(f+2),0,0.1) ,tr); Row[j]=1; Row=paste(Row,collapse=' ')
      u = c(u,Row)}
    para = c(para,u)
  }
  u = 
'OPTION SNP_file geno
OPTION use_yams
OPTION EM-REML 10
OPTION approx_loglike
OPTION maxrounds 20'
  para = c(para,u)
  write.table(para,'par',sep='',row.names=F,col.names=F,quote=F)
  # Part 5
  cat('done\n')
}
compiler::cmpfun(ram)

read_blup = function(){
  uu = 'trait effect level  solution'
  x = read.table('solutions',sep = 'k')
  levels(x$V1) = c(levels(x$V1),uu)
  x[1,1] = 'trait effect level  solution'
  write.table(x,'solutions2',row.names=F,col.names=F,quote = F)
  x = read.table('solutions2',T)
  xx = read.table('name_index')
  q = x$solution[x$effect==max(x$effect)]
  q = matrix(q,ncol = max(x$trait),byrow = TRUE)
  # for(i in 1:max(x$trait)) q[,i]=q[,i]+x$solution[i]
  rownames(q) = xx$V1
  write.table(q,'solutions2', col.names=F,quote = F)
  return(q)
}
compiler::cmpfun(read_blup)

######################

# Run first example?
if(FALSE){
  # Get example data
  data(soynam,package='SoyNAM')
  data = data.line; names(data)[4] = 'ID'
  form = yield | R8 ~ environ + ID
  gen = gen.raw[1:500,]; # 500 individuals for this example
  data = droplevels(data[which(data$family<5),])
  rm(data.check,data.line,gen.raw)
  # Generate files
  ram(form,data,gen)
  # RUN
  system('airemlf90.exe',input='par')
  # Collect results
  blups = read_blup() 
  plot(blups)
}

######################

# Run second example?
if(FALSE){
  # Get example data
  require(NAM)
  data(met)
  # Add spatial covariates
  Obs = data.frame(Obs,spY=round(SPC(Obs$YLD,Obs$Block,Obs$Row,Obs$Block),2))
  Obs = data.frame(Obs,spACC=round(SPC(Obs$ACC,Obs$Block,Obs$Row,Obs$Block),2))
  Obs = data.frame(Obs,spDTM=round(SPC(Obs$DTM,Obs$Block,Obs$Row,Obs$Block),2))
  Gen = data.matrix(Gen[1:200,])
  Obs = droplevels.data.frame(Obs[which(Obs$ID%in%rownames(Gen)),])
  # Generate files
  ram(form=YLD|ACC|DTM~Year+ID+spY+spACC+spDTM,data=Obs,gen=Gen)
  # Run
  system('airemlf90.exe',input='par')
  # Collect results
  blups = data.frame(read_blup())
  plot(blups)
}


