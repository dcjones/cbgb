
library(hexbin)


a <- read.table( 'intron_exon_read_density.tab', header=FALSE )
names(a) <- c( 'id', 'exon_reads', 'intron_reads', 'exon_len', 'intron_len' )

b <- subset( a, exon_reads/exon_len < .5 & intron_reads/intron_len < .5 &
             intron_len > 1000 & exon_len > 100 )

c <- subset( b, (exon_reads/exon_len) > 0.1 & (intron_reads/intron_len) < 0.2 )
#attach(c)
#reg <- glm( intron_reads ~ intron_len + I(exon_reads/exon_len),
            #family=poisson('identity'), start = c(0,1,0.1) )
#summary(reg)
#warnings()
#detach(c)



reg <- lm( I(c$exon_reads/c$exon_len) ~ I(c$intron_reads/c$intron_len) - 1 )
reg_slope <- reg[[1]][[1]]
print(reg)
#print(reg[[1]][[1]])
#print(reg[[1]][[2]])




# simulated data
s <- rnbinom( length(b$V2), mu=0.1*(b$exon_reads/b$exon_len), size=b$intron_len )


colramp <- function(n) { return(plinrain(n,beg=2)) }


png( 'exon_intron_read_density.png', width=800, height=800 )



hexbinplot(
     I(exon_reads/exon_len) ~ (intron_reads/intron_len),
     data=b,
     trans=log,
     inv=exp,
     xbins=250,
     colramp=colramp,
     aspect=1,
     panel=function(...){
         panel.fill('black')
         panel.hexbinplot(...)


         # standard deviations
         reg_stddev <- function( c, x ) {
             mu <- x*reg_slope
             return( mu + c*sqrt(mu)/10e10 )
         }
         ys <- 0.001 * 0:1000

         xs <- mapply( reg_stddev, 1, ys )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #xs <- mapply( reg_stddev, 1, ys )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #ys <- mapply( reg_stddev, 2, xs )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #ys <- mapply( reg_stddev, 3, xs )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #ys <- mapply( reg_stddev, -1, xs )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #ys <- mapply( reg_stddev, -2, xs )
         #panel.xyplot( xs, ys, type='l', col='white' )

         #panel.abline( reg, col='white' )
     },
     xlab='intronic',
     ylab='exonic'
     )


dev.off()
