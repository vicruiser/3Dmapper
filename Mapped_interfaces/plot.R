library(ggplot2)
# Bar plot
plot_table <- fread("/home/vruizser/Escritorio/PDB_CALC_DIST_finalRes/mpi_v89_5_vs_v94.csv")

tb <- ddply(plot_table, .(mapped.pi),summarize,
                        chain = sum(chain),
                        chain.1 = sum(chain.1),
                        ensembl.prot.id = sum(ensembl.prot.id),
                        sstart = sum(sstart),
                        send = sum(send),
                        real.pos = sum(real.pos)
            )

(counts <- with(diamonds, table(cut, clarity)))
#            clarity
# cut           I1  SI2  SI1  VS2  VS1 VVS2 VVS1   IF
# Fair         210  466  408  261  170   69   17    9
# Good          96 1081 1560  978  648  286  186   71
# Very Good     84 2100 3240 2591 1775 1235  789  268
# Premium      205 2949 3575 3357 1989  870  616  230
# Ideal        146 2598 4282 5071 3589 2606 2047 1212

(counts <- with(plot_table, table(mapped.pi, chain)))

ggplot(diamonds, aes(clarity, fill = cut)) + 
  geom_bar(position = 'identity', alpha = .3)