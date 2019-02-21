t1=c(1,2,3,4,5,2,1)
u=unique(t1)
u2=u*5
inx=sapply(t1, function(i) which(i==u) )

t2=u2[unlist(inx)]



gene=  "ENSG00000093072"

rp2 <- mclapply(1:nrow(rp.r), function(i) t(rbind(rp.f,rp.r[i,]))) ## get ref panel fsnps and rsnp in the same format I have functions from simulations
names(rp) <- rs.full$id

urp=unique(rp2)

## calculate P(H|G) ref panel for each rsnp, split to avoid memory issues 
mod= length(urp) %% block
block=400
int=as.integer(length(urp)/block)

urp.hap.pairs=list()

i=0
while(i*block< length(urp)){
  if(i*block+mod == length(urp)) {
    #urp.hap.pairs <- c(urp.hap.pairs, i*block +1, (i+1)*block+mod)
    urp.hap.pairs[[(i+1)]] = mclapply(urp[(i*block+1):((i+1)*block+mod)], p.hap.pair)
  } else {
    urp.hap.pairs[[(i+1)]]  = mclapply(urp[(i*block+1):((i+1)*block)], p.hap.pair)
  #urp.hap.pairs <- c(urp.hap.pairs, i*block+1, (i+1)*block)
  }
 
  i=i+1
  print(i)
  gc()
}
  

## memory issue, very big object which I only need temp, sorted by processing downstream steps in the same call.