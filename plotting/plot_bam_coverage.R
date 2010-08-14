
# A script to produce a plots of the broad coverage of a set of bam files across
# all chromosomes.
#
# Daniel Jones <dcjones@cs.washington.edu>
# May 11, 2010
#


# USE:
# From the command line,
# $ R --no-save < plot_bam_coverage.R
# (this may take a several minutes or more)


# PREREQUISITES
# The following packages are required:
# 1. GenomicRanges (BAM file manipulation) 
# 2. Rsamtools     (BAM file manipulation)
# 2. gtools        (sorting by chromosome name)
# 3. Cairo         (pretty graphics)
# 4. ggplot2       (pretty graphics)
#
# 
# 1,2 requires at least R version 2.11 and at least Bioconductor 2.6.
# They can be installed with the following command.
#  > source('http://bioconductor.org/biocLite.R')
#  > biocLite( c('Rsamtools','GenomicRanges') )
#
# The last three you can get from CRAN: 
#  > install.packages( c('gtools','Cairo','ggplot2') )


# These are the BAM files that will be plotted. Change this variable as
# appropriate.
sets <- Sys.glob('*/hits.bam')

# So the plot is readable, reads are binned into regions of size m.
# A smaller m makes a more detailed plot, but possibly with lines to fine to
# see without a high powered microscope.
m <- 100000



# Part I: Count reads in each bin from each BAM file.
library(GenomicRanges) 

b <- list()


for( set in sets ) {
    a <- readGappedAlignments( set )
    a <- narrow( a, end=-width(a) )
    c <- coverage(a)

    for( chr in names(c) )
    {
        v <- as.integer(c[[chr]])
        for( i in 1:(ceiling(length(v)/m)) ) {
            while( i > length(b[[chr]]) ) {
                b[[chr]] <- c( b[[chr]], 0 )
            }

            v_sums <- sum( v[((i-1)*m+1):(i*m)], na.rm=TRUE )
            b[[chr]][[i]] <- b[[chr]][[i]] + v_sums
        }
    }
}


# Part II: Rearrange data into a data frame.
d <- data.frame()
d['seq'] <- NULL
d['pos'] <- NULL
d['count'] <- NULL


for( chr in names(c) ) {
    tmp <- as.data.frame( b[[chr]] )
    colnames(tmp) <- c('count')
    tmp['pos'] <- m * as.integer(rownames(tmp))
    tmp['seq'] <- chr
    d <- rbind(d,tmp)
}

library(gtools)

d <- d[ mixedorder(d$seq), ]
d$seq <- factor(d$seq, levels=unique(d$seq))


# Part III: Make a pretty plot.
library(Cairo)
library(ggplot2)

CairoPNG('coverage.png', width=1100, height=2000)


p <- qplot( pos, (count), data=d, geom='area' )
p <- p + scale_x_continuous(name='Position')
p <- p + scale_y_continuous(name='Read Count')
p <- p + opts( panel.grid.major = theme_blank(),
               panel.grid.minor = theme_line(size=0.1,colour='white') )
p <- p + facet_grid( seq ~ . )

print(p)


dev.off()


