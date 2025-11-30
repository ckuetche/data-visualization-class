#----plotting a volcano plot for question 3 in PS3.----#

#WARNING: make sure to download data from galaxy as .tsv, not tabular (r cannot read it)

#set the working directory
setwd("/storage/group/exd44/pracbio25/students/rfk5434/problemset3/question3")

#Load libraries, data, and verify-----------------------------------------------

#load your required library
library(ggplot2)

#read the data by specifying the file name, header, separator, and stringsAsFactors (false)
results <- read.delim("results_file_attempt2.tsv", header=FALSE, stringsAsFactors=FALSE)

#add the headers manually using colnames(), refer to galaxy
colnames(results) <- c("GeneID", "baseMean", "log2FoldChange", "StdErr", "Wald-Stats", "P-value", "P-adj")

head(results) #verify 34642 obs of 7 variables

#the y axis is -log10(p-value), so add that as a column - use backticks for literal column names
results$negLogPadj <- -log10(results$`P-adj`)
head(results) #verify 34642 obs of 8 variables

#Make the base scatter plot-----------------------------------------------------

#log2FoldChange on the x axis, and -log10(p-value) for the y axis with aes
ggplot(data=results, aes(x=log2FoldChange, y=negLogPadj)) + #+ adds layers to the plot
  geom_point() #defines how data is represented as points

#the next modifications will use this structure, so save it as a variable
vol_plot <- ggplot(data=results, aes(x=log2FoldChange, y=negLogPadj)) +
  geom_point()

vol_plot

#plot quandrants----------------------------------------------------------------

#assume significant genes have a p_value<0.05 & 0.5 > log2FC > 2
#y intercept will isolate genes with p_value<0.05, and the x intercept will filter downregulated (log2FC 0.5)
#and upregulated genes (log2FC 2) genes. use geom_hline or geom_vline (hor, ver) and specify linetype (dashed)

vol_plot <- vol_plot +
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +
  geom_vline(xintercept=c(log2(0.5),log2(2)),
             linetype="dashed")
vol_plot
#update your variable at the end

#adding color, size, and transparency-------------------------------------------
#label genes as upregulated and downregulated (fc>=2, fc<0.5, both (p<0.05))
#use dplyr::case_when, which checks conditions in a vector for each element simultaneously
library(dplyr)
#create a new datastruct to categorize by regulation
modified_results <- results %>%
    mutate(
        regulation=case_when(
            log2FoldChange >= 2 & `P-adj` < 0.05 ~ "Upregulated",
            log2FoldChange <= 0.5 & `P-adj` < 0.05 ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )
#mutate() will take your data, create or modify new columns, and return the dataset, %>% is a pipe_operator

#obtain gene counts
modified_results %>% 
    count(regulation)
#from modified results, keep the first occurrence of the regulation column, extract that single column as a vector
#check your gene categories
modified_results %>%
  distinct(regulation) %>%
  pull()
cols <- c("Upregulated"="#E64B35", "Downregulated"="#4DBBDB", "Not Significant"="grey")
sizes <- c("Upregulated"=2, "Downregulated"=2, "Not Significant"=1)
alphas <- c("Upregulated"=1, "Downregulated"=1, "Not significant"=0.5)
#alphas - 1 is fully opaque, 0 is fully transparent

modified_results %>%
  ggplot(aes(x=log2FoldChange,
              y=negLogPadj,
              fill=regulation,
              size=regulation,
              alpha=regulation)) +
  geom_point(shape=21,
             color="black") + #outline color
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +
  geom_vline(xintercept=c(log2(0.5),log2(2)),
             linetype="dashed") +
  scale_fill_manual(values=cols) + #modify point color
  scale_size_manual(values=sizes) + #modify point size
  scale_alpha_manual(values=alphas) + ##scale_*_manual() are manual override functions we use
  guides(
    size="none", #hide redundant legend
    alpha="none"
  ) +
  theme_minimal() +
  labs( #add titles - labs stands for labels
    title="Volcano Plot of Differential Expression",
    x="Log2(Fold Change)",
    y="-Log10(Adjusted P-value)",
    fill="regulation" #set legend title
  )



