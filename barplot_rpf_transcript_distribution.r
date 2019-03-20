#!/usr/bin/env Rscript

require(data.table)
require(ggplot2)
in1="./data/output/rpf_5p_density/miseq_data.rpf_in_regs.miseq_data"
in2="./data/output/rpf_5p_density/miseq_data.rpf_in_regs.miseq_data_norm"
output_dir="./data/output/figures"

mydt1 <- fread(in1)
mydt2 <- fread(in2)

mydt1[,Type:='Counts']
mydt2[,Type:='RPK']

mydt <- rbind(mydt1,mydt2)

mydt[,Sample:= gsub('run_umi_(.*)_umied',"\\1",Sample)]
mysum <- melt(mydt, id.vars=c('Sample','Type'))
myplot <- ggplot(mysum, aes(x=Sample, y=value, fill=variable))
myplot <- myplot + geom_bar(stat='identity', position='stack',width=0.8) + theme_bw() +
    ylab('Average % of reads') + xlab(NULL) + facet_wrap(~Type) 
myplot
ggsave(paste(outdir, 'rpf_in_regs.pdf', sep=""),myplot)
