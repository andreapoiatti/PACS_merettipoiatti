points.2D.data =function(SpacePoints,Data,zlim = NA,Nx=100,Ny=100, samePlot = FALSE, nlevels = 100){
  # colore = rev(rainbow(nlevels,start=0,end=4/6))
  colore = heat.colors(nlevels)
  spessore = seq(from = 1, to = 3, length= nlevels)
  
  if(is.na(zlim)[1]){
    Min = min(Data)
    Max = max(Data)
  }else{
    Min = zlim[1]
    Max = zlim[2]
  }
  
  if(samePlot)
  {
    GenMatrix = matrix(data = Data, nrow = 1, ncol = dim(SpacePoints)[1] ,byrow=T)
    Gen = GenMatrix[1,]
    
    for(pippo in 1:nlevels){
      points(SpacePoints[which(round((Gen-Min)/(Max-Min)*nlevels)==pippo),1],SpacePoints[which(round((Gen-Min)/(Max-Min)*nlevels)==pippo),2],
             col=colore[pippo],cex=spessore[pippo],pch=16)
    }
  }
  else
  {
    xmin = min(SpacePoints[,1]) - abs(min(SpacePoints[,1]))/20
    xmax = max(SpacePoints[,1]) + abs(max(SpacePoints[,1]))/20
    
    ymin = min(SpacePoints[,2]) - abs(min(SpacePoints[,2]))/20
    ymax = max(SpacePoints[,2]) + abs(max(SpacePoints[,2]))/20
    
    xGen = seq(xmin,xmax,length.out = Nx)
    yGen = seq(ymin,ymax,length.out = Ny)
    
    
    
    PosMatrix = expand.grid(xGen,yGen)
    PlotVector = rep(0, dim(PosMatrix)[1])
    PlotMatrix = matrix(data = PlotVector, nrow = length(xGen), ncol = length(yGen),byrow=F)
    
    GenMatrix = matrix(data = Data, nrow = 1, ncol = dim(SpacePoints)[1] ,byrow=T)
    Gen = GenMatrix[1,]
    
    image(xGen,yGen,PlotMatrix,col="white",xlab="",ylab="", zlim = c(Min,Max),asp = 1)
    
    for(pippo in 1:nlevels){
      points(SpacePoints[which(round((Gen-Min)/(Max-Min)*nlevels)==pippo),1],SpacePoints[which(round((Gen-Min)/(Max-Min)*nlevels)==pippo),2],
             col=colore[pippo],cex=spessore[pippo],pch=16)
    }
  }
  
}