#!/bin/Rscript

library(tidyverse)
library(cowplot)
library(grid)

options(scipen = 999)
args <- commandArgs(TRUE)


# read custom lib
{
  bin_path_ <- strsplit(args[1],"=")[[1]][2]
  source(paste0(bin_path_,"/custom_R_func_lib_.R"))
}

# read args
{
  opts_ <- parse_args(args)
  
  accession_ <- opts_$accession_ %||% stop("accession_ required")
  query_file_ <- opts_$query_file_ %||% stop("query_file_ required")
  soloLTR_file_ <- opts_$soloLTR_file_ %||% stop("soloLTR_file_ required")
  karyo_file_ <- opts_$karyo_file_ %||% stop("karyo_file_ required")
}

print(accession_)

# get input data
{
  # get query data
  file_i_ <- query_file_
  print(query_file_)
  if (file.info(file_i_)$size > 0) {
    query_data_ <- read.table(query_file_,header = F,sep = "\t",skip = 1,comment.char = "")
    colnames(query_data_) <- c("queryID_","FL_length_","lLTR_length_","rLTR_length_","intdom_length_","motif_","length_diff_","tg_ca_","cd_hit_")
    query_data_ <- query_data_ %>%
      mutate(class_ = gsub(".*#(.*?)/.*","\\1",queryID_))
  }else{
    print(paste0("Missing ", accession_,"_wga_LTRRT_lengths.txt file!"))
    quit(status = 1)
  }
  
  # get soloLTR data
  file_i_ <- soloLTR_file_
  print(soloLTR_file_)
  if (file.info(file_i_)$size > 0) {
    soloLTR_data_ <- read.table(soloLTR_file_,header = F,sep = "\t",skip = 1,comment.char = "")
    colnames(soloLTR_data_) <- c("seqID_","start_","end_","queryID_","strand_","piden_","locus_","locus_type_","TSD_alg_","BLAST_single_hit_","homology_","TSD_")
  }else{
    print(paste0("Missing ", accession_,"_wga_LTR_BLAST_overlap.txt file!"))
    quit(status = 1)
  }
  
  # get karyo data
  file_i_ <- karyo_file_
  print(karyo_file_)
  if (file.info(file_i_)$size > 0) {
    karyo_data_ <- read.table(karyo_file_,header = F,sep = "\t")
    colnames(karyo_data_) <- c("seqID_","seqID_length_")
  }else{
    print(paste0("Missing ", accession_,"_karyotype_ext_ file!"))
    quit(status = 1)
  }
  
}

# stats
{  
  gs_ <- karyo_data_$seqID_length_ %>% sum()
  
  total_LTRRT_n_ <- query_data_ %>% nrow
  total_LTRRT_space_ <- query_data_$FL_length_ %>% sum
  total_LTR_lib_ <- query_data_[which(query_data_$length_diff_ == "PASS" & query_data_$tg_ca_ == "PASS" & query_data_$cd_hit_ == "PASS"),] %>% nrow
  total_LTR_lib_prop_ <- (total_LTR_lib_/total_LTRRT_n_)*100 
  total_LTR_lib_prop_ <- paste0(round(total_LTR_lib_prop_,digits = 1),"%")
  total_solo_ltr_n_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>% nrow
  total_solo_ltr_bp <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
    mutate(length_ = (end_-start_)+1,) %>%
    select(length_) %>%
    sum()
  
  main_stats_ <- data.frame("values_" = c(gs_,total_LTRRT_n_,total_LTRRT_space_,total_LTR_lib_,total_LTR_lib_prop_,total_solo_ltr_n_,total_solo_ltr_bp))
  main_stats_$labels_  <- c("genome size (bp)","total LTRRT","LTRRT space (bp)","LTR library size", "compression","total soloLTR","soloLTR space (bp)")
  main_stats_$labels_ <- main_stats_$labels_ %>%
    {factor(.,levels = unique(main_stats_$labels_))}
  
  p_stats_ <- ggplot() +
    geom_point(data = main_stats_,aes(x = 0, y = labels_), color = NA) +
    geom_label(data = main_stats_[grep("soloLTR",main_stats_$labels_),],aes(x = 0, y = labels_, label = values_),hjust = 0,border.color = NA,label.r = unit(0, "lines"),fill = "#FFF8DC", size.unit = "pt", size = 10) +
    geom_label(data = main_stats_[-grep("soloLTR",main_stats_$labels_),],aes(x = 0, y = labels_, label = values_),hjust = 0,border.color = NA,label.r = unit(0, "lines"),fill = NA, size.unit = "pt", size = 10) +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_discrete(limits = rev) +    
    coord_cartesian(xlim = c(-0.05,1.05),clip = "off") +
    labs(title = accession_) +
    theme_void() +
    theme(axis.text.y = element_text(hjust = 1),
          plot.margin = margin(r =10,l = 10,t = 2.5,b = 2.5,unit = "pt"),
          plot.title = element_text(face = "bold", hjust = 0, size = 12),
          text = element_text(size = 10, color = "black"))
}  

# query length dist
# soloLTR length dist
{
  length_dist_FL_df_ <- query_data_ %>%
    select(queryID_,FL_length_) %>%
    distinct() %>%
    mutate(feature_ = "Full Length LTR-RT",
           length_ = FL_length_,
           ID_ = queryID_) %>%
    select(!c(FL_length_,queryID_)) %>%
    mutate(mean_ = mean(length_,na.rm = T),
           median_ = median(length_,na.rm = T))
  
  length_dist_lLTR_df_ <- query_data_ %>%
    select(queryID_,lLTR_length_) %>%
    distinct() %>%
    mutate(feature_ = "long terminal repeat",
           length_ = lLTR_length_,
           ID_ = queryID_) %>%
    select(!c(lLTR_length_,queryID_))  %>%
    mutate(mean_ = mean(length_,na.rm = T),
           median_ = median(length_,na.rm = T))
  
  length_dist_soloLTR_df_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
    mutate(feature_ = "soloLTR",
           length_ = (end_-start_)+1,
           ID_ = "queryID_") %>%
    select(feature_,length_,ID_) %>%
    mutate(mean_ = mean(length_,na.rm = T),
           median_ = median(length_,na.rm = T))
  
  length_dist_df_ <- rbind(length_dist_FL_df_,length_dist_lLTR_df_,length_dist_soloLTR_df_)
  
  maxVal_ <- max(length_dist_df_$length_[which(length_dist_df_$feature_ != "Full Length LTR-RT")],na.rm = T)
  binwidth_ <- 10
  p_dist_ <- ggplot(length_dist_df_) +
    annotate(geom = "point", y = 0, x = maxVal_, color = NA) +
    geom_vline(aes(xintercept = mean_), linetype = "solid") +
    geom_vline(aes(xintercept = median_), linetype = "dotted") +
    geom_histogram(data = length_dist_df_, aes(x = length_,fill = feature_),binwidth = binwidth_, color = NA) +
    xlab("Length (bp)") +
    ylab("Frequency") +
    scale_fill_manual(values = c("Full Length LTR-RT" = "darkblue","long terminal repeat" = "darkblue","soloLTR" = "orangered")) +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(expand = expansion(.025)) +
    facet_wrap(vars(feature_), ncol = 1, scales = "free") +
    coord_cartesian(clip = "off", expand = T) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.ontop = T,
          axis.title = element_text(),
          strip.clip = "off",
          legend.position = "none", 
          text = element_text(color = "black"))
}

# TSD library diversity
{
  tsd_lib_sum_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
    group_by(TSD_alg_) %>%
    summarise(n_ =n()) %>%
    group_by(n_) %>%
    ungroup() %>%
    mutate(bin_n_ = case_when(n_ == 1 ~ 1,
                              n_ == 2 ~ 2,
                              n_ == 3 ~ 3,
                              n_ <= 10 & n_ > 3 ~ 4,
                              n_ > 10 ~ 5),
           bin_n_ = as.character(bin_n_)) %>%
    group_by(bin_n_) %>%
    summarise(freq_unique_ = n(),
              freq_all_ = sum(n_)) %>%
    mutate(unique_ = freq_unique_/sum(freq_unique_),
           all_ = freq_all_/sum(freq_all_))
    
  tsd_lib_sum_anno_unique_ <- tsd_lib_sum_ %>%
    mutate(cumsum_ = cumsum(unique_)) %>%
    mutate(midpoint_ = ifelse(bin_n_ == "1", cumsum_/2, (unique_/2)+dplyr::lag(cumsum_)))
    
  tsd_lib_sum_anno_all_ <- tsd_lib_sum_ %>%
    mutate(cumsum_ = cumsum(all_)) %>%
    mutate(midpoint_ = ifelse(bin_n_ == "1", cumsum_/2, (all_/2)+dplyr::lag(cumsum_)))
  
  palette_ <- data.frame(hex_color_ = c("#FF4500","#FF7A41","#FFA674","#FFD0A7","#FFF8DC"))
  palette_$breaks_ <- c("1","2","3","4","5")
  palette_$labels_ <- c("1","2","3","[4,10]",">10")
  
  p_tsd_div_ <- ggplot() +
    geom_col(data = tsd_lib_sum_, aes(x = "0_unique", y = unique_, fill = bin_n_), color = "white",linewidth = 1, position = position_stack(reverse = T), width = .815) +
    geom_col(data = tsd_lib_sum_, aes(x = "1_all", y = all_, fill = bin_n_), color = "white",linewidth = 1, position = position_stack(reverse = T), width = .815) +
    geom_text(data = tsd_lib_sum_anno_unique_, aes(x = "0_unique", y = midpoint_, label = freq_unique_), size.unit = "pt", size = 6.5, color = "black") +
    geom_text(data = tsd_lib_sum_anno_all_, aes(x = "1_all", y = midpoint_, label = freq_all_), size.unit = "pt", size = 6.5, color = "black") +
    scale_fill_manual(values = palette_$hex_color_, breaks = palette_$breaks_, labels = palette_$labels_) +
    scale_x_discrete(expand = expansion(mult = c(.5, 0)),labels = ~ gsub(".*_","",.)) +
    scale_y_continuous(breaks = seq(0,1,0.25), labels = ~ ifelse(. != 0, paste0(.*100,"%"), 0),expand = expansion(0)) +
    coord_cartesian(clip = "off") +
    labs(title = "TSD diversity") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          plot.title.position = "plot",
          plot.title = element_text(color = "black", face = "bold",hjust = 0),
          legend.title = element_blank(),
          legend.position = "right",
          legend.direction = "vertical",
          legend.justification = "center",
          text = element_text(color = "black"))

  # TSD frequency
    tsd_lib_freq_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
      group_by(TSD_alg_) %>%
      summarise(n_ =n()) %>%
      group_by(n_) %>%
      ungroup() %>%
      mutate(bin_n_ = case_when(n_ == 1 ~ 1,
                                n_ == 2 ~ 2,
                                n_ == 3 ~ 3,
                                n_ <= 10 & n_ > 3 ~ 4,
                                n_ > 10 ~ 5),
             bin_n_ = as.character(bin_n_))
    
    p_tsd_freq_ <- ggplot() +
      annotate(geom = "rect", xmin = 0.25,xmax = 1.75, ymin = 1, ymax = 2, color = NA, fill = palette_$hex_color_[1]) +
      annotate(geom = "rect", xmin = 0.25,xmax = 1.75, ymin = 2, ymax = 3, color = NA, fill = palette_$hex_color_[2]) +
      annotate(geom = "rect", xmin = 0.25,xmax = 1.75, ymin = 3, ymax = 4, color = NA, fill = palette_$hex_color_[3]) +
      annotate(geom = "rect", xmin = 0.25,xmax = 1.75, ymin = 4, ymax = 11, color = NA, fill = palette_$hex_color_[4]) +
      annotate(geom = "rect", xmin = 0.25,xmax = 1.75, ymin = 11, ymax = Inf, color = NA, fill = palette_$hex_color_[5]) +
      geom_boxplot(data = tsd_lib_freq_, aes(x = 1, y = n_),width = .65,fill = "white",outliers = T) +
      scale_y_log10(expand = expansion(0),
                    breaks = 10^(0:5),
                    minor_breaks = as.vector(
                      outer(10^(0:4), 2:9)),
                    guide = guide_axis(minor.ticks = TRUE)) +
      scale_color_manual(values = palette_$hex_color_, breaks = palette_$breaks_, labels = palette_$labels_) +
      labs(title = " ") +
      ylab("Frequency (Log10)") +
      coord_cartesian(xlim = c(-0.1,2.1),clip = "off") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_line(),
            axis.ticks.y = element_line(),
            axis.ticks.length.y = unit(5, "pt"),
            axis.minor.ticks.length.y = unit(2, "pt"),
            plot.title.position = "plot",
            plot.title = element_text(color = "black", face = "bold",hjust = 1),
            legend.title = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.justification = "center")
}

# mirroring comparison between intact and solo
# visual intact to solo ratio
# soloLTR
soloLTR_sum_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
  mutate(feature_ = "soloLTR",
         length_ = (end_-start_)+1,
         class_ = gsub(".*#(.*?)/.*","\\1",queryID_)) %>%
  group_by(class_,feature_) %>%
  summarise(n_ = n(),
            space_ = sum(length_)) %>%
  mutate(n_prop_ = n_/sum(n_),
         space_prop_ = space_/sum(space_))

# all features together
soloLTR_sum_all_ <- soloLTR_data_[which(soloLTR_data_$TSD_ == "PASS"),] %>%
  mutate(feature_ = "soloLTR",
         length_ = (end_-start_)+1,
         class_ = "all") %>%
  group_by(class_,feature_) %>%
  summarise(n_ = n(),
            space_ = sum(length_)) %>%
  mutate(n_prop_ = n_/sum(n_),
         space_prop_ = space_/sum(space_))

soloLTR_sum_ <- rbind(soloLTR_sum_,soloLTR_sum_all_)

# LTRRT
LTRRT_sum_ <- query_data_ %>%
  mutate(feature_ = "Full Length LTR-RT",
         length_ = FL_length_,
         class_ = gsub(".*#(.*?)/.*","\\1",queryID_)) %>%
  group_by(class_,feature_) %>%
  summarise(n_ = n(),
            space_ = sum(length_)) %>%
  mutate(n_prop_ = n_/sum(n_),
         space_prop_ = space_/sum(space_))

# all features together
LTRRT_sum_all_ <- query_data_ %>%
  mutate(feature_ = "Full Length LTR-RT",
         length_ = FL_length_,
         class_ = "all") %>%
  group_by(class_,feature_) %>%
  summarise(n_ = n(),
            space_ = sum(length_)) %>%
  mutate(n_prop_ = n_/sum(n_),
         space_prop_ = space_/sum(space_))

LTRRT_sum_ <- rbind(LTRRT_sum_,LTRRT_sum_all_)


# LTR
lLTR_sum_ <- query_data_ %>%
  mutate(feature_ = "long terminal repeat",
         length_ = lLTR_length_,
         class_ = gsub(".*#(.*?)/.*","\\1",queryID_)) %>%
  group_by(class_,feature_) %>%
  summarise(n_ = n(),
            space_ = sum(length_)) %>%
  mutate(n_prop_ = n_/sum(n_),
         space_prop_ = space_/sum(space_))


# n_
{
  step_ <- 0.20
  # soloLTR
  n_start_ <- soloLTR_sum_$n_[which(soloLTR_sum_$class_ == "all")] %>% 
    {sum(.)*seq(0,1,step_)}
  n_end_ <- soloLTR_sum_$n_[which(soloLTR_sum_$class_ == "all")] %>% 
    {c(sum(.)*seq(step_,1,step_),1)}
  breaks_ <- c(seq(step_,1,step_),1)
  label_ <- paste0(breaks_*100,"%")
  soloLTR_prop_ <- cbind(n_start_,n_end_,breaks_,label_) %>%
    data.frame()
  soloLTR_prop_[,c(1,2,3)] <- apply(soloLTR_prop_[,c(1,2,3)],2,as.numeric)
  soloLTR_prop_ <- soloLTR_prop_[-nrow(soloLTR_prop_),]
  
  # LTRRT
  n_start_ <- LTRRT_sum_$n_[which(LTRRT_sum_$class_ == "all")] %>% 
    {sum(.)*seq(0,1,step_)}
  n_end_ <- LTRRT_sum_$n_[which(LTRRT_sum_$class_ == "all")] %>% 
    {c(sum(.)*seq(step_,1,step_),1)}
  breaks_ <- c(seq(step_,1,step_),1)
  label_ <- paste0(breaks_*100,"%")
  LTRRT_prop_ <- cbind(n_start_,n_end_,breaks_,label_) %>%
    data.frame()
  LTRRT_prop_[,c(1,2,3)] <- apply(LTRRT_prop_[,c(1,2,3)],2,as.numeric)
  LTRRT_prop_ <- LTRRT_prop_[-nrow(LTRRT_prop_),]
  
  # ratio
  ratio_ <- rbind(LTRRT_sum_,soloLTR_sum_) %>%
    group_by(class_) %>%
    summarise(ratio_ = n_[feature_ != "soloLTR"]/n_[feature_ == "soloLTR"])
  
  
  grey_palette_ <- rep(c("grey94","grey98"),nrow(LTRRT_prop_))
  names(grey_palette_) <- LTRRT_prop_$label_
  
  maxVal_l_ <- max(LTRRT_prop_$n_end_,na.rm = T)*-1
  maxVal_r_ <- max(soloLTR_prop_$n_end_,na.rm = T)
  prop_label_ <- max(soloLTR_prop_$n_start_,na.rm = T)+((max(soloLTR_prop_$n_end_,na.rm = T)-max(soloLTR_prop_$n_start_,na.rm = T))/2)
 
  step_ <- 0.20
  p_ratio_ <- ggplot() +
    geom_rect(data = LTRRT_prop_[,], aes(xmin = n_start_*-1, xmax = n_end_*-1, ymin = -Inf, ymax = Inf, fill = label_)) +
    geom_rect(data = soloLTR_prop_[,], aes(xmin = n_start_, xmax = n_end_, ymin = -Inf, ymax = Inf, fill = label_)) +
    annotate(geom = "text",x = prop_label_, y = Inf, label = paste0("\n\n",step_*100,"%"), size.unit = "pt", size = 6,color = "grey45") +
    geom_col(data = LTRRT_sum_, aes(x = n_*-1, y = class_), color = NA, fill = "darkblue", width = .605) +
    geom_col(data = soloLTR_sum_, aes(x = n_, y = class_), color = NA, fill = "orangered", width = .605) +
    geom_text(data = ratio_,aes(y = class_, x = 0+(maxVal_r_*0.115), label = round(ratio_,digits = 2)), size.unit = "pt", size = 4.85,hjust = 1) +
    geom_vline(xintercept = 0) +
    xlab("Obs (#)") +
    facet_wrap(vars(class_),ncol = 1, scales = "free") +
    scale_x_continuous(labels = ~ ifelse(. < 0, .*-1,.),expand = expansion(0)) +
    scale_fill_manual(values = grey_palette_) +
    coord_cartesian(clip = "on",xlim = c(maxVal_l_,maxVal_r_)) +
    theme(panel.background = element_rect(color = NA, fill = NA),
          panel.grid = element_blank(),
          panel.ontop = T,
          axis.title.y = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          strip.clip = "off",
          strip.background = element_rect(color = NA, fill = NA),
          strip.text = element_text(hjust = 0), 
          text = element_text(color = "black"))
  
  }

# karyotyping
{
  
  # only sequences with at least one soloLTR hit will be printed
  seqID_i_list_ <- soloLTR_data_$seqID_[which(soloLTR_data_$TSD_ == "PASS")] %>% unique()
  seqID_i_list_ <- karyo_data_$seqID_[which(karyo_data_$seqID_ %in% seqID_i_list_)] 

  maxSeqID_ <- max(karyo_data_$seqID_length_[which(karyo_data_$seqID_ %in% seqID_i_list_)],na.rm = T)
  linewidth_ <- .25
  p_LTRRT_list_ <- list()

  for(i_ in 1:length(seqID_i_list_)){
    seqID_i_ <- seqID_i_list_[i_]
    seqID_length_i_ <- karyo_data_$seqID_length_[which(karyo_data_$seqID_ == seqID_i_)]
    print(seqID_i_)
    # intact sum
    if(nrow(query_data_[grep(paste0(seqID_i_,"_"),query_data_$queryID_),]) == 0){
      next
    }
    binwidth_ <- 100000
    LTRRT_dist_sum_ <- query_data_[grep(paste0(seqID_i_,"_"),query_data_$queryID_),] %>%
      mutate(seqID_ = seqID_i_, 
             start_ = as.numeric(gsub(".*_([0-9]+)_[0-9]+#.*/.*","\\1",queryID_)),
             end_ = as.numeric(gsub(".*_[0-9]+_([0-9]+)#.*/.*","\\1",queryID_)),
             class_ = gsub(".*#(.*?)/.*","\\1",queryID_)) %>%
      mutate(dist_start_ = findInterval(start_,seq(1,seqID_length_i_,length.out = round(seqID_length_i_/binwidth_)),all.inside = T),
             length_ = (end_-start_)+1) %>%
      group_by(dist_start_) %>%
      summarise(n_ = n(),
                space_ = sum(length_),
                space_prop_ = space_/binwidth_)
    
    # soloLTR
    if(nrow(soloLTR_data_[which(soloLTR_data_$seqID_ == seqID_i_ & soloLTR_data_$TSD_ == "PASS"),]) == 0){
      next
    }
    soloLTR_loc_ <- soloLTR_data_[which(soloLTR_data_$seqID_ == seqID_i_ & soloLTR_data_$TSD_ == "PASS"),] %>%
      mutate(class_ = gsub(".*#(.*?)/.*","\\1",queryID_)) %>%
      select(seqID_,start_,end_,class_,TSD_alg_) %>%
      distinct() %>%
      mutate(dist_start_ = findInterval(start_,seq(1,seqID_length_i_,length.out = round(seqID_length_i_/binwidth_)),all.inside = T),
             length_ = (end_-start_)+1) %>%
      group_by(dist_start_,class_) %>%
      summarise(n_ = n(),
                space_ = sum(length_),
                space_prop_ = space_/binwidth_)
    
    max_occupancy_ <- range(soloLTR_loc_$space_prop_)[2]
    
    step_ <- 1
    dist_start_ <- seq(1,seqID_length_i_/binwidth_,step_) %>%
      data.frame()
    colnames(dist_start_) <- "dist_start_"
    LTRRT_dist_sum_ <-  merge(LTRRT_dist_sum_,dist_start_,all.y = T)
    LTRRT_dist_sum_[,c(2,3,4)] <- apply(LTRRT_dist_sum_[,c(2,3,4)],2,function(x) ifelse(is.na(x) == T, 0,x))
    
    class_list_ <- unique(query_data_$class_)
    
    p_LTRRT_ <- ggplot() +
      annotate(geom = "rect", xmin = 0, xmax = seqID_length_i_, ymin = 0, ymax = 1, fill = "grey96") +
      geom_hline(yintercept = 1,color = "white") +
      geom_hline(yintercept = .50,color = "white") +
      geom_area(data = LTRRT_dist_sum_, aes(x = dist_start_*binwidth_, y = space_prop_), linewidth = linewidth_, fill = "darkblue", color = NA, alpha = .2) +
      geom_line(data = LTRRT_dist_sum_, aes(x = dist_start_*binwidth_, y = space_prop_), linewidth = linewidth_, color = "darkblue")
    
    gap_ <- 0.4
    for(j_ in 1:length(class_list_)){
      class_j_ <- class_list_[j_]
      y_ <- (seq(gap_,10,gap_)[j_])*-1
      soloLTR_loc_class_ <- soloLTR_loc_[which(soloLTR_loc_$class_ == class_j_),]
      soloLTR_loc_class_$y_ <- y_
      p_LTRRT_ <- p_LTRRT_ +
        annotate(geom = "rect", xmin = 0, xmax = seqID_length_i_, ymin = y_+0.025, ymax = y_-0.025, fill = "grey96", color = NA) +
        geom_tile(data = soloLTR_loc_class_, aes(x = dist_start_*binwidth_, y = y_, fill = space_prop_), color = NA, height = .35, width = binwidth_)
    }
    
    y_ <- (seq(gap_,10,gap_)[j_+1])*-1
    y_ <- y_+(gap_*0.25)
    
    n_breaks_ <- 4
    breaks_x_ <- scales::breaks_extended(n = n_breaks_)(c(0, seqID_length_i_)) %>%
      {.[. <= seqID_length_i_]} %>%
      data.frame()
    colnames(breaks_x_)[c(1)] <- "value_"
    
    p_LTRRT_ <- p_LTRRT_ +
      annotate(geom = "segment", x = 0, xend = seqID_length_i_, y = y_, yend = y_, color = "black") +
      geom_segment(data = breaks_x_, aes(x = value_, xend = value_, y = y_, yend = y_-0.165)) +
      scale_y_continuous(breaks = c((seq(gap_,10,gap_)[1:j_])*-1,0,.5,1), labels = c(paste0(rev(class_list_)," "),"0 ","50% ","100% "),expand = expansion(mult = c(0,.025)))  +
      scale_x_continuous(expand = expansion(0), breaks = breaks_x_$value_, label = breaks_x_$value_) +
      scale_fill_gradient(low = "#FFA674", high = "orangered") +
      labs(title = seqID_i_,
           subtitle = paste0("max bin val: ",round(max_occupancy_*100),"%")) +
      coord_cartesian(xlim = c(0,maxSeqID_),clip = "off") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_text(size = 10,color = "black"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 8,color = "black"),
            axis.line.x = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(t=5,r=5,b=0,l=5,"pt"),
            plot.title = element_text(hjust = 0,face = "bold"),
            plot.subtitle = element_text(color = "black",size = 7.5,hjust = 0),
            text = element_text(color = "black", size = 10))
    
    p_LTRRT_list_[[length(p_LTRRT_list_)+1]] <- p_LTRRT_
  }
  
}


# grid
{
  top_right_ <- cowplot::plot_grid(p_dist_,p_tsd_freq_,p_tsd_div_,ncol = 3,rel_widths = c(1,0.425,.775), align = "hv", axis = "lb",labels = c("B","D"))  
    
  p_ratio_<- cowplot::plot_grid(NULL,p_ratio_,ncol = 2,rel_widths = c(.05,1))
  top_left_ <- cowplot::plot_grid(p_stats_,p_ratio_,ncol = 1,rel_heights = c(0.525,1),labels = c("A","C"))
    
  top_ <- cowplot::plot_grid(top_left_,top_right_,ncol = 2,rel_widths = c(.485,1),align = "h", axis = "b")
    
  p_karyo_ <- cowplot::plot_grid(plotlist=p_LTRRT_list_,ncol = 1,align = "hv",axis = "tbrl")
  p_karyo_ <- cowplot::plot_grid(p_karyo_,NULL, rel_heights = c(1,0.0001), ncol = 1, labels = c("E",""))
  
}

# print
{
  n_karyo_ <- length(seqID_i_list_)
  karyo_r_ <- (n_karyo_*15)/9.25
  top_r_ <- 5.35
  
  pdf(paste0(accession_,"_soloLTRseeeker_summmary_.pdf"), width = 10.5, height = top_r_+karyo_r_)
  grid.newpage()
  pushViewport(viewport(
    layout = grid.layout(
      2, 1,
      heights = unit.c(unit(top_r_, "in"), unit(karyo_r_, "in"))
    )
  ))
  pushViewport(viewport(layout.pos.row = 1))
  print(top_, vp = viewport())
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2))
  print(p_karyo_, vp = viewport())
  popViewport(2)
  dev.off()
}
