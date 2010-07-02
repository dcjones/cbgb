
library(Cairo)
library(ggplot2)


plot_hashed_reads <- function( in_fn, out_fn ) {

    a <- read.table( in_fn, header=F )
    a <- a[order(a$V1,decreasing=T),]
    n <- length(a[,1])
    a <- cbind(a,seq(n))


    names(a) <- c('count','seq','i')

    #a$count <- a$count / sum(a$count)

    CairoPNG( out_fn, width=1000, height=800 )

    p <- ggplot( a, aes( x=i, y=cumsum(count) ) )
    p <- p + geom_line()
    print(p)

    dev.off()
}
