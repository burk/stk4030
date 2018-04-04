display.matrix = function(x,nrow=100,ncol=100, ...)
{
    ##
    ## display a matrix as an image with correct layout and autoscaling
    ## read dimentions from x if x is a matrix, otherwise use nrow and ncol
    ##
    
    mi = min(x)
    ma = max(x)

    if (is.matrix(x) && all(dim(x)>1))
    {
        n = dim(x)
    }
    else
    {
        x = matrix(as.vector(x),nrow=nrow,ncol=ncol)
        n = dim(x)
    }
    y = x
    for(j in 1:n[2]) y[1:n[1],j] = x[n[1]:1,j]
    image((t(y)-mi)/(ma-mi), col=gray(seq(0,1,len=256), ...))
}

