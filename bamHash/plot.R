
library(Cairo)
library(ggplot2)

CairoFonts(
    regular="Univers LT Std,Univers LT Std 55:style=55 Light,Regular",
    bold="Univers LT Std,Univers LT Std 45 Light:style=45 Roman,Regular"
    )


plot_hashed_reads <- function( in_fn, out_fn ) {

    a <- read.table( in_fn, header=F )
    a <- a[order(a$V1,decreasing=T),]
    n <- length(a[,1])
    a <- cbind(a,seq(n))


    names(a) <- c('seq','count','i')

    #a$count <- a$count / sum(a$count)

    CairoPNG( out_fn, width=1000, height=800 )

    p <- ggplot( a, aes( x=i, y=cumsum(count) ) )
    p <- p + geom_line()
    print(p)

    dev.off()
}
