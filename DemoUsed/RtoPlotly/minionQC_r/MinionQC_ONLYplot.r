#!/usr/bin/env Rscript

# MinIONQC version 1.3.5
# Copyright (C) 2017 onwards Robert Lanfear
#
# For license see https://github.com/roblanf/minion_qc/blob/master/LICENSE

#only plot, only 1D QC, for RmarkdownTest

# supress warnings
options(warn=-1)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plotly))


#p1m = 1.0
p1m = 0.5

# this is how we label the reads at least as good as q
q_title = paste("Q>=", q, sep="")

# build the map for R9.5
p1 = data.frame(channel=33:64, row=rep(1:4, each=8), col=rep(1:8, 4))
p2 = data.frame(channel=481:512, row=rep(5:8, each=8), col=rep(1:8, 4))
p3 = data.frame(channel=417:448, row=rep(9:12, each=8), col=rep(1:8, 4))
p4 = data.frame(channel=353:384, row=rep(13:16, each=8), col=rep(1:8, 4))
p5 = data.frame(channel=289:320, row=rep(17:20, each=8), col=rep(1:8, 4))
p6 = data.frame(channel=225:256, row=rep(21:24, each=8), col=rep(1:8, 4))
p7 = data.frame(channel=161:192, row=rep(25:28, each=8), col=rep(1:8, 4))
p8 = data.frame(channel=97:128, row=rep(29:32, each=8), col=rep(1:8, 4))

q1 = data.frame(channel=1:32, row=rep(1:4, each=8), col=rep(16:9, 4))
q2 = data.frame(channel=449:480, row=rep(5:8, each=8), col=rep(16:9, 4))
q3 = data.frame(channel=385:416, row=rep(9:12, each=8), col=rep(16:9, 4))
q4 = data.frame(channel=321:352, row=rep(13:16, each=8), col=rep(16:9, 4))
q5 = data.frame(channel=257:288, row=rep(17:20, each=8), col=rep(16:9, 4))
q6 = data.frame(channel=193:224, row=rep(21:24, each=8), col=rep(16:9, 4))
q7 = data.frame(channel=129:160, row=rep(25:28, each=8), col=rep(16:9, 4))
q8 = data.frame(channel=65:96, row=rep(29:32, each=8), col=rep(16:9, 4))

map = rbind(p1, p2, p3, p4, p5, p6, p7, p8, q1, q2, q3, q4, q5, q6, q7, q8)

add_cols <- function(d, min.q){
    # take a sequencing sumamry file (d), and a minimum Q value you are interested in (min.q)
    # return the same data frame with the following columns added
        # cumulative.bases and cumulative.bases.time
        # hour of run
        # reads.per.hour

    d = subset(d, mean_qscore_template >= min.q)

    if(nrow(d)==0){
        flog.error(paste("There are no reads with a mean Q score higher than your cutoff of ", min.q, ". Please choose a lower cutoff and try again.", sep = ""))
        quit()
    }

    d = merge(d, map, by="channel")

    d = d[with(d, order(as.numeric(start_time))), ] # sort by start time
    d$cumulative.bases.time = cumsum(as.numeric(d$sequence_length_template))

    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    d$hour = d$start_time %/% 3600

    # add the reads generated for each hour
    reads.per.hour = as.data.frame(table(d$hour))
    names(reads.per.hour) = c("hour", "reads_per_hour")
    reads.per.hour$hour = as.numeric(as.character(reads.per.hour$hour))
    d = merge(d, reads.per.hour, by = c("hour"))
    return(d)
}


load_summary <- function(filepath, min.q){
    # load a sequencing summary and add some info
    # min.q is a vector of length 2 defining 2 levels of min.q to have
    # by default the lowest value is -Inf, i.e. includes all reads. The
    # other value in min.q is set by the user at the command line
    d = read_tsv(filepath, col_types = cols_only(channel = 'i',
                                                num_events_template = 'i',
                                                sequence_length_template = 'i',
                                                mean_qscore_template = 'n',
                                                sequence_length_2d = 'i',
                                                mean_qscore_2d = 'n',
                                                start_time = 'n',
                                                calibration_strand_genome_template = 'c'))

    # remove the control sequence from directRNA runs
    if("calibration_strand_genome_template" %in% names(d)){
        d = subset(d, calibration_strand_genome_template != "YHR174W")
    }

    if("sequence_length_2d" %in% names(d)){
        # it's a 1D2 or 2D run
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_2d))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_2d))
        d$num_events_template = NA
        d$start_time = as.numeric(as.character(d$start_time))

    }else{
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_template))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_template))
        d$num_events_template = as.numeric(as.character(d$num_events_template))
        d$start_time = as.numeric(as.character(d$start_time))
    }

    d$events_per_base = d$num_events_template/d$sequence_length_template

    flowcell = basename(dirname(filepath))

    # add columns for all the reads
    d1 = add_cols(d, min.q[1])
    d1$Q_cutoff = "All reads"

    # add columns for just the reads that pass the user Q threshold
    d2 = add_cols(d, min.q[2])
    d2$Q_cutoff = q_title

    # bind those two together into one data frame
    d = as.data.frame(rbindlist(list(d1, d2)))

    # name the flowcell (useful for analyses with >1 flowcell)
    d$flowcell = flowcell

    # make sure this is a factor
    d$Q_cutoff = as.factor(d$Q_cutoff)

    keep = c("hour","start_time", "channel", "sequence_length_template", "mean_qscore_template", "row", "col", "cumulative.bases", "cumulative.bases.time", "reads_per_hour", "Q_cutoff", "flowcell", "events_per_base")
    d = d[keep]

    return(d)
}


log10_minor_break = function (...){
    # function to add minor breaks to a log10 graph
    # hat-tip: https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        minor_breaks =
            rep(log10(seq(1, 9, by=1)), times = n_major)+
            rep(major_breaks, each = 9)
        return(10^(minor_breaks))
    }
}

log10_major_break = function (...){
    # function to add major breaks to a log10 graph
    # hat-tip: https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        return(10^(major_breaks))
    }
}

binSearch <- function(min, max, df, t = 100000) {
    # binary search algorithm, thanks to https://stackoverflow.com/questions/46292438/optimising-a-calculation-on-every-cumulative-subset-of-a-vector-in-r/46303384#46303384
    # the aim is to return the number of reads in a dataset (df)
    # that comprise the largest subset of reads with an N50 of t
    # we use this to calculte the number of 'ultra long' reads
    # which are defined as those with N50 > 100KB
    mid = floor(mean(c(min, max)))
    if (mid == min) {
        if (df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[min]/2))] < t) {
            return(min - 1)
        } else {
            return(max - 1)
        }
    }

    n = df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[mid]/2))]
    if (n >= t) {
        return(binSearch(mid, max, df))
    } else {
        return(binSearch(min, mid, df))
    }
}


channel.summary <- function(d){
    # calculate summaries of what happened in each of the channels
    # of a flowcell
    a = ddply(d, .(channel),
              summarize,
              total.bases = sum(sequence_length_template),
              total.reads = sum(which(sequence_length_template>=0)),
              mean.read.length = mean(sequence_length_template),
              median.read.length = median(sequence_length_template))
    b = melt(a, id.vars = c("channel"))
    return(b)
}


single.flowcell <- function(input.file, output.dir, q=7, base.dir = NA){
    # wrapper function to analyse data from a single flowcell
    # input.file is a sequencing_summary.txt file from a 1D run
    # output.dir is the output directory into which to write results
    # q is the cutoff used for Q values, set by the user
    # base.dir is the base directory if and only if the user supplied a base directory
    # we use base.dir to name flowcells in a sensible way
    flog.info(paste("Loading input file:", input.file))
    d = load_summary(input.file, min.q=c(-Inf, q))

    flowcell = unique(d$flowcell)

    # output goes with the sequencing summary file unless otherwise specified
    
        	
    d$muxes = seq(from = 0, to = max(d$hour), by = 8)

    # set up variable sizes
    if(smallfig == TRUE){ p1m = 0.5 }else{ p1m = 1.0 }
    if(smallfig == TRUE){ p2m = 0.6 }else{ p2m = 1.0 }


    return(d)
}


#main()
d = single.flowcell(input.file, output.dir, q)
p1 = ggplot(d, aes(x = sequence_length_template, fill = Q_cutoff)) +
        geom_histogram(bins = 300) +
        scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") +
        theme(text = element_text(size = 15)) +
        xlab("Read length (bases)") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)

p2 = ggplot(d, aes(x = mean_qscore_template, fill = Q_cutoff)) +
        geom_histogram(bins = 300) +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") +
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)

p3 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x=start_time/3600, y=sequence_length_template, colour = mean_qscore_template)) +
        geom_point(size=1.5, alpha=0.35) +
        scale_colour_viridis() +
        labs(colour='Q')  +
        scale_y_log10() +
        facet_grid(row~col) +
        theme(panel.spacing = unit(0.5, "lines")) +
        xlab("Hours into run") +
        ylab("Read length")# +
        #theme(text = element_text(size = 40), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), legend.text=element_text(size=18), legend.title=element_text(size=24))
p5 = ggplot(d, aes(x=start_time/3600, y=cumulative.bases.time, colour = Q_cutoff)) +
        geom_vline(xintercept = d$muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        geom_line(size = 1) +
        xlab("Hours into run") +
        ylab("Total yield in bases") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
p6 = ggplot(d, aes(x=sequence_length_template, y=cumulative.bases, colour = Q_cutoff)) +
        geom_line(size = 1) +
        xlab("Minimum read length (bases)") +
        ylab("Total yield in bases") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
    xmax = max(d$sequence_length_template[which(d$cumulative.bases > 0.01 * max(d$cumulative.bases))])
p6 = p6 + scale_x_continuous(limits = c(0, xmax))

p7 = ggplot(d, aes(x=start_time/3600, y=sequence_length_template, colour = Q_cutoff, group = Q_cutoff)) +
        geom_vline(xintercept = d$muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() +
        xlab("Hours into run") +
        ylab("Mean read length (bases)") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        ylim(0, NA)

p8 = ggplot(d, aes(x=start_time/3600, y=mean_qscore_template, colour = Q_cutoff, group = Q_cutoff)) +
        geom_vline(xintercept = d$muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() +
        xlab("Hours into run") +
        ylab("Mean Q score") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        ylim(0, NA)        

f = d[c("hour", "reads_per_hour", "Q_cutoff")]
f = f[!duplicated(f),]
g = subset(f, Q_cutoff=="All reads")
h = subset(f, Q_cutoff==q_title)
max = max(f$hour)
# all of this is just to fill in hours with no reads recorded
all = 0:max
add.g = all[which(all %in% g$hour == FALSE)]
if(length(add.g)>0){
        add.g = data.frame(hour = add.g, reads_per_hour = 0, Q_cutoff = "All reads")
        g = rbind(g, add.g)
    }
add.h = all[which(all %in% h$hour == FALSE)]
if(length(add.h)>0){
        add.h = data.frame(hour = add.h, reads_per_hour = 0, Q_cutoff = q_title)
        h = rbind(h, add.h)
    }
i = rbind(g, h)
i$Q_cutoff = as.character(i$Q_cutoff)
i$Q_cutoff[which(i$Q_cutoff==q_title)] = paste("Q>=", q, sep="")
p9 = ggplot(i, aes(x=hour, y=reads_per_hour, colour = Q_cutoff, group = Q_cutoff)) +
        geom_vline(xintercept = d$muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_point() +
        geom_line() +
        xlab("Hours into run") +
        ylab("Number of reads per hour") +
        ylim(0, NA) +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads"))        

p10 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) +
        geom_point(alpha=0.05, size = 0.4) +
        scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) +
        labs(colour='Events per base\n(log scale)\n')  +
        theme(text = element_text(size = 15)) +
        xlab("Read length (bases)") +
        ylab("Mean Q score of read")
if(max(d$events_per_base, na.rm=T)>0){
        # a catch for 1D2 runs which don't have events per base
        p10 = p10 + scale_colour_viridis(trans = "log", labels = scientific, option = 'inferno')
    }
# we keep it a bit wider, because the legend takes up a fair bit of the plot space

c = channel.summary(subset(d, Q_cutoff=="All reads"))
c10 = channel.summary(subset(d, Q_cutoff==q_title))
c$Q_cutoff = "All reads"
c10$Q_cutoff = q_title
cc = rbind(c, c10)
cc$variable = as.character(cc$variable)
cc$variable[which(cc$variable=="total.bases")] = "Number of bases per channel"
cc$variable[which(cc$variable=="total.reads")] = "Number of reads per channel"
cc$variable[which(cc$variable=="mean.read.length")] = "Mean read length per channel"
cc$variable[which(cc$variable=="median.read.length")] = "Median read length per channel"
p11 = ggplot(cc, aes(x = value, fill = Q_cutoff)) + geom_histogram(bins = 30) +
        facet_grid(Q_cutoff~variable, scales = "free_x") +
        #theme(text = element_text(size = 15), axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(text = element_text(size = 8), axis.text.x = element_text(angle = 60, hjust = 1)) +
        guides(fill=FALSE) +
        scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +
        guides(fill=FALSE)