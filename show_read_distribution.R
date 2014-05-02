# Overview of new preprocessing
library(reshape2)
library(ggplot2)
library(scales)

# Input values for every sample as seperated vector
# structure c(RAW, paired trimmed, single trimmed, concat, no_dup) 
sample60 <- c(RAW = 1401439, p_trim = 1072089, s_trim = 217047, 
              concat = 490160, no_dup =  1870846)
sample62 <- c(RAW = 6155862, p_trim = 4696630, s_trim = 883440, 
              concat = 1760627, no_dup = 8514557)
sample64 <- c(RAW = 8743868, p_trim = 3173163, s_trim = 7271064, 
              concat = 2612057, no_dup = 5550285)
sample66 <- c(RAW = 17654788, p_trim = 5236025, s_trim = 5996006, 
              concat = 3816133, no_dup = 12641540)
sample68 <- c(RAW = 1736051, p_trim = 1366893, s_trim = 246074, 
              concat = 564401, no_dup = 2415265)
sample70 <- c(RAW = 1991046, p_trim = 1469638, s_trim = 311022, 
              concat = 511392, no_dup = 2738325)
sample74 <- c(RAW = 11769965, p_trim = 3434548, s_trim = 1805583, 
              concat = 2858577, no_dup = 5813667)
sample76 <- c(RAW = 1682298, p_trim = 1191996, s_trim = 284656, 
              concat = 400040, no_dup = 2268302)
sample78 <- c(RAW = 1713223, p_trim = 1216594, s_trim = 272004, 
              concat = 310785, no_dup = 2388324)
sample80 <- c(RAW = 1585176, p_trim = 1154371, s_trim = 247431, 
              concat = 398013, no_dup = 2157281)
sample82 <- c(RAW = 1233661, p_trim = 874699, s_trim = 194189, 
              concat = 333192, no_dup = 1609940)

# combine in one vector
data <- data.frame(sample60, sample62, sample64, 
                   sample66, sample68, sample70, 
                   sample74, sample76, sample78, 
                   sample80, sample82)
# rearrenge vector for ggplot2
data2 <- melt(data)
# new factor for legend 
desc <- factor(rep(c("RAW", "paired end trimmed", "single end trimmed",
                     "concatinated", "single without duplicates"),
                   length=nrow(data2)),
              levels=c("RAW", "paired end trimmed", "single end trimmed",
                       "concatinated", "single without duplicates"))
# new vector for facet_wrap
dest <- factor(c(rep("free living", 25),rep("aqua culture", 30)),
               levels=c("free living", "aqua culture"))
# combine to new data.frame
df <- cbind(data2, desc, dest)
# build plot
pdf("graphs/read_distribution_meta-pipeline.pdf")
  ggplot(df, aes(x = variable, y = value, fill = desc)) + 
    # type of plot
    geom_bar(stat = "identity", position = "dodge") +
    # names of the axis
    xlab("\nSample") + ylab("Number of Reads")  +
    # seperate into two windows
    theme_bw() + facet_wrap(~dest, scales="free_x") + 
    # adjust Y-Axis
    scale_y_continuous(labels = comma, breaks = pretty_breaks(n=15)) +
    # change colors
    scale_fill_manual(values = c("brown2", "chartreuse4", "chartreuse1", 
                                "darkorchid1", "#6495ED")) + 
    # change legend title
    guides(fill = guide_legend("Type of reads")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Distribution of Reads after preprocessing steps\n")
dev.off()
