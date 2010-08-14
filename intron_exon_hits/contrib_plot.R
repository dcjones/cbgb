
library(lattice)

a <- read.table( 'intron_exon_read_density.disjoint.tab', header=FALSE )
names(a) <- c( 'id', 'exon_reads', 'intron_reads', 'exon_len', 'intron_len' )




sum_reads <- a$exon_reads + a$intron_reads
#sum_reads <- a$intron_reads

sum_reads <- rev(sum_reads[order(sum_reads)])
sum_reads <- sum_reads[ sum_reads > 0 ]

total <- sum(sum_reads)

print(head(sum_reads))

xs <- 1:length(sum_reads)
ys <- mapply( function(n) { return(sum(sum_reads[1:n])) }, xs )

png( 'exonintron_contrib.png', width=800, height=800 )
xyplot( I(ys/total) ~ I(xs/length(sum_reads)), type='l',
            panel = function(...) {
                panel.xyplot(...)
                panel.grid(h=-10,v=-10)
            } )

dev.off()
