#library(e1071)
library(rpart)
#library(caTools)
#library(preprocessCore)
library(stats)
#library(rJava)
#library(xlsxjars)
library(openxlsx)
library(ggplot2)
library(vioplot)
library(dendextend)
library(stringr)
library(pvclust)
library(stats)
#library(DESeq)
#library(DESeq2)
#library(oligo)
#library(CONOR)
#library(HARMONY)
libarary(dplyr)
library(magrittr)

#rm(list=ls(all=TRUE))
#memory.size (max = TRUE)
#memory.limit ( size = 4000)

bar_y_scale = function(dendrogram) {
    return(attr(dendrogram, "height") * 0.1)
}

###
## criteria for evaluating the clustering quality ###
###

# output:
# list with the plot and significant area
evaluate_hclust <-
    function(
        tree,
        # a hclust object
        classes,
        # a character vector of class labels, order should be the same as the data
        random_trials = 300,
        # how much random trials (label shuffling OR random tree generation) to perform
        seed = 451,
        # random seed
        progress = TRUE
        ) {
        
    require(easypackages)
    ## will install all the required packages if they aren't present
    suppressMessages(packages('magrittr', 'dplyr', 'stringr', 'dendroextras', 'progress', prompt = TRUE))
    # these limits are used both in dendro_ig and in this function
    possible_clusts <- nrow(tree$merge) + 1
    cut_limits <- c(2, possible_clusts - 1) 

    ## creating a progress bar for convenience
    if (progress) {
    pb <-
        progress_bar$new(
            total = random_trials,
            format = "Running random trials... [:bar] :percent Time left::eta"
            )
    }
    ## this function calculates information gain trajectory from a single tree
    dendro_ig <- function(tree, classes, progress = FALSE) {
    # options.ig <- list(lowerp = 0, upperp = 1) # deprecated
    stopifnot(inherits(tree, 'hclust'))
    ref_ext <- 0
    
    ## function which performs metric calculation
    infgain <- function(real_classes, progress = progress) {
        ## helper function to calculate ig for a given cut
        ig_cut <- function(cut_classes, real_classes, progress = progress) {
            ## entropies and relative frequencies of all the clusters
            cluster_entropies <-
                numeric(n_distinct(cut_classes)) %>%
                set_names(unique(cut_classes))
            
            cluster_prevalences <-
                table(cut_classes) %>%
                divide_by(sum(.))
            
            ref <- 0
            ## another helper; calculates single cluster's entropy
            entropy_in_cluster <- function(clust_classes) {
                ## an obvious case
                if (length(clust_classes) == 1) return(0)
                
                clust_classes %>%
                    table() %>%
                    ## there shouldn't be zeros in the table
                    divide_by(sum(.)) %>%
                    inset(. == 0, 1) %>%
                    ## entropy is measured in bits
                    {-1*.*log2(.)} %>%
                    sum()
            }
            for (clust in unique(cut_classes)) {
                cluster_entropies[as.character(clust)] <- entropy_in_cluster(real_classes[cut_classes == clust])
                # print(clust)
                # print(cluster_entropies[as.character(clust)])
            }
            cluster_norm_entr <- -cluster_entropies * cluster_prevalences
            # H0 = 
            # ref <- entropy_in_cluster(real_classes)
            # variable in original function's envirorment for normalizing info_gains
            # the upper one is "more local"
            if (ref_ext == 0) ref_ext <<- entropy_in_cluster(real_classes)
            # print(str_c('Ref:', ref))
            # print(str_c('Normalized entropies: ', str_c(cluster_norm_entr, collapse = ', ')))
            ## calculating information gain
            IG <-
                cluster_norm_entr %>%
                sum() %>%
                add(ref_ext)
            # naming accordingly
            IG
        }
        ## a vector with information gains as by [all specified] tree cut levels
        ## the first element will remain zero
        info_gains <- numeric(cut_limits[2])
        ## iterating along [all specified] cuts, 2 to N-1
        for (cut_depth in cut_limits[1]:cut_limits[2]) {
            # the old way, should be slower, left for historical reasons
            # info_gains[cut_depth] <- ig_cut(tree, cutree(tree, k = cut_depth), real_classes)
            # it should find all_cuts_tree matrix two environments up
            info_gains[cut_depth] <- ig_cut(all_cuts_tree[cut_depth, ], real_classes)
            # an unfinished prototype, may speed things up just a bit more
            # info_gains <- sapply(all_cuts_tree, ig_cut, real_classes = real_classes)
        }
        ## adding start (above) and end points, they are always strictly 0 and 1
        # normalized information gain
        info_gains %<>% 
            divide_by(ref_ext) %>%
            c(., 1)
        # this is just for the random trial progress to show up
        if (progress) pb$tick()
        return(info_gains)
    }
    
    out <-
        infgain(classes, progress = progress) %>%
        ## we should know for which clusters the output was obtained
        set_names(c(1, cut_limits[1]:cut_limits[2], possible_clusts))
    out
    }
    
    possible_clusts <- nrow(tree$merge) + 1
    actual_clusters <- 1:possible_clusts
    # calculating tree cuts just one time, a matrix is supposedly faster than a list
    # each row represents a cut, from top to bottom
    all_cuts_tree <- matrix(0, nrow = possible_clusts, ncol = possible_clusts)
    for (cut_depth in cut_limits[1]:cut_limits[2]) {
        all_cuts_tree[cut_depth, ] <- cutree(tree, k = cut_depth)
    }
    ## data for the true labels
    igs <- dendro_ig(tree = tree, classes = classes)
    
    ## list of shuffled class labels
    set.seed(seed)
    rtrials <- replicate(random_trials, sample(classes), simplify = FALSE)
    ## a matrix of information gains, each column represents a random trial
    ## each run returns a numeric vector (2...n_samples - 1) of total length n_samples - 2
    # random trials, vapply is for safety + a bit of speed
    trial_igs <-
        vapply(
            rtrials,
            dendro_ig,
            tree = tree,
            progress = progress,
            FUN.VALUE = numeric(length(classes))
            )
    
    # packing the result
    data_rtrials <-
        data_frame(
            k = rep(actual_clusters, random_trials),
            igain = as.numeric(trial_igs),
            run = rep(1:random_trials, each = possible_clusts)
        )

    # IG(1) = 0, IG(possible_clusts) = 1
    # cut_limits <- c(2, possible_clusts - 1)
    ## writing results to data frames in order to later plot them and calculate various statistics
    data_igs <-
        data_frame(
            k = actual_clusters,
            igain = igs
        )
    ## a list with two data frames, for original and shuffled data
    list_out <-
        list(
            # observed data
            data_igs = data_igs,
            # random data
            data_rtrials = data_rtrials,
            # useful for calculating ideal splits
            class_labels = classes
            # determines fot color on the watermelon plot later
        )
    return(list_out)
}

# helper function used both for area calculation and for plotting the ideal curve
ideal_ig_trajectory <- function(class_labels) {
    # entropy for a single cluster
    entropy <- function(labels) {
        ## an obvious case
        if (length(labels) == 1) return(0)
        
        labels %>%
            table() %>%
            ## there shouldn't be zeros in the table
            divide_by(sum(.)) %>%
            inset(. == 0, 1) %>%
            ## entropy is measured in bits
            {-1*.*log2(.)} %>%
            sum()
    }
    # order in which the classes will be excluded; ties don't matter
    class_freqs <-
        class_labels %>%
        table() %>%
        sort(decreasing = TRUE) %>%
        # frequencies
        divide_by(sum(.))

    # base entropy
    H0 <- entropy(class_labels)
    # a vector containing ideally-splitted classes' entropies
    max_ig_trajectory <-
        c(0, rep(1, length(class_labels) - 1)) %>%
        set_names(
            c(
                # root
                "root",
                # removed class
                names(class_freqs),
                # additional stuff
                rep("", length(class_labels) - length(class_freqs) - 1)
            )
        )
    levels_excluded <- character(0)

    for (level in names(class_freqs)[1:(length(class_freqs) - 1)]) {
        levels_excluded %<>% c(level)
        max_ig_trajectory[level] <- 1 - class_freqs[level]*entropy(class_labels[!class_labels %in% levels_excluded])
    }
    # 1 if everything has been grouped correctly
    names(max_ig_trajectory)[1:length(class_freqs) + 1] %<>% str_c("removed ", .)
    return(max_ig_trajectory)
}

## a separate function to obtain the criterion value
## takes an output of evaluate_hclust
watermelon_area <- function(ig_datalist) {
    data_rtrials <- ig_datalist$data_rtrials
    data_igs <- ig_datalist$data_igs
    data_classes <- ig_datalist$class_labels
    cut_limits <- c(1, max(data_rtrials$k))
    
    rquantiles <-
        data_rtrials %>%
        group_by(k) %>%
        dplyr::summarise(igain = quantile(igain, 0.95)) %$%
        igain

    ideal_traj <- ideal_ig_trajectory(data_classes)
    ## histogram-like estimate of total area,
    ## negative area also counts!
    area <- sum(data_igs$igain - rquantiles) / sum(ideal_traj - rquantiles)
    return(area)
}

## a separate function to plot the result with more customization
## takes an output of evaluateHclust
watermelon_plot <-
    function(
        ig_datalist,
        # an output of evaluateHclust
        plot_name = 'Watermelon plot',
        # plot name displayed in the heading
        bootstrap_ci = TRUE,
        # include bootstrap 95%CI? will take some time to calculate
        bootstrap_runs = 500,
        # how much bootstrap trials should be performed?
        dots_percent = 0.6,
        # percent of dots shown on the plot for faster plotting
        seed = 451,
        # random seed
        show_legend = TRUE
        ) {
    require(easypackages)
    suppressMessages(libraries('tidyr', 'dplyr', 'forcats', 'ggplot2', 'magrittr', 'stringr'))
    
    data_igs <-
        ig_datalist$data_igs %>%
        mutate(`Curve type` = 'Observed data')
    data_rtrials <- ig_datalist$data_rtrials
    
    data_ideal <-
        data_frame(
            k = 1:nrow(data_igs),
            igain = ideal_ig_trajectory(ig_datalist$class_labels),
            `Curve type` = 'Ideal class separation'
        )

    data_criterion <-
        data_rtrials %>%
        group_by(k) %>%
        dplyr::summarise(igain = quantile(igain, 0.95)) %>%
        mutate(`Curve type` = '95% of random trials') %>%
        bind_rows(., data_igs, data_ideal) %>%
        mutate(`Curve type` = factor(`Curve type`))
    
    ## used to draw the area
    ## transformed to satisfy row-wise ymin / ymax easthetics format
    data_area <-
        data_criterion %>%
        filter(`Curve type` %in% c("Observed data", "95% of random trials")) %>%
        tidyr::spread(`Curve type`, igain) %>%
        mutate(area_sign = ifelse(`Observed data` > `95% of random trials` | `95% of random trials` == 1, '+', '-') %>% factor() %>% fct_rev()) # %>%
        # if this isn't done this will result in a geom_ribbon stretching glitch
        # filter(`95% of random trials` != `Observed data`)
    
    ## for faster plotting
    if (dots_percent != 1) {
        ## randomly sampling the dots
        ## the plot will look more chaotic, the outliers appear quite variable at each k
        set.seed(seed)
        data_rtrials %<>% filter(run %in% 1:floor(max(run)*dots_percent))
    }
    
    area <- watermelon_area(ig_datalist)
    
    ## dot size should be dependent on the number of samples (the graph becomes smaller)
    ## a heuristic
    minor_dot_size <- 1.4 * (4.5/sqrt(nrow(data_igs)))
    major_dot_size <- 2.65 * (4.5/sqrt(nrow(data_igs)))
    # a stringie which contains confidence interval
    boot_string <- ci_string <- character(0)
    if (bootstrap_ci) {
        # reuqires bootstrap.R to be sourced first!
        wm_bootstrap <-
            watermelon_bootstrap(
                ig_datalist = ig_datalist, 
                step = -1, 
                bootstrap_trials = bootstrap_runs,
                seed = seed
                )
        wm_ci_df <- watermelon_area_ci(wm_bootstrap)
        wm_ci <- round(c(wm_ci_df$lower_bound, wm_ci_df$upper_bound), 4)
        
        # pieces of information that'd be plotted
        ci_string <- str_c('\n95% bootstrap CI: [', wm_ci[1], '; ', wm_ci[2], ']')
        boot_string <- str_c('\nNumber of bootstrap resamples = ', bootstrap_runs)
    }

    main_colors <-
        c(
            'Ideal class separation' = 'dodgerblue3',
            'Observed data' = 'forestgreen',
            '95% of random trials' = 'firebrick2'
        )

    plot_ig <-
        ggplot(data_igs, aes(x = k, y = igain)) +
        geom_abline(
            slope = 1 / (max(data_igs$k) - 1),
            intercept = -1 / (max(data_igs$k) - 1),
            color = 'darkgrey',
            size = 2.5
        ) +
        geom_jitter(
            data = data_rtrials,
            mapping = aes(x = k, y = igain),
            color = 'firebrick2',
            alpha = 0.07,
            size = major_dot_size,
            width = 0.15,
            height = 0
        ) +
        geom_line(
            data = data_criterion,
            aes(x = k, y = igain, color = `Curve type`),
            size = 2.0,
            inherit.aes = FALSE
        ) +
        scale_color_manual(values = main_colors) + 
        geom_point(shape = 16, size = minor_dot_size, color = 'black') +
        theme(
            legend.position = c(0.8, 0.2),
            legend.background = element_rect(linetype = 'solid', color = 'black', fill = 'lightgrey'),
            panel.background = element_blank(), 
            panel.grid.major = element_line(color = 'lightgrey'),
            panel.grid.minor = element_line(color = 'lightgrey')
                ) + 
        # guides(color = guide_legend(override.aes = list(shape = c(5, 5, 5)), reverse = FALSE)) + 
        labs(x = 'Dendrogram depth',
                 y = 'Normalized Information Gain',
                 title = str_c(
                     plot_name,
                     '\n',
                     str_c('Statistically significant area = ',
                                 round(area, 4),
                                 ' units; ',
                                 ci_string)
                     ),
                 caption = str_c(
                     'Number of randomizations = ', max(data_rtrials$run),
                     # either empty or contains info on a number of bootstrap trials conducted
                     boot_string,
                     '\nSeed = ', seed
                     )
                 )
                if(!show_legend) plot_ig <- plot_ig + guides(color = FALSE)
    
    plot_ig
}


# writes the results in separate txt files
write_wm <- function(wm, path) {
    field_names <- c("data_igs", "data_rtrials", "class_labels")
    for (field_name in names(wm)) {
        write_tsv(tbl_df(wm[[field_name]]), str_c(path, field_name, ".txt"))
    }
}

# reads the results from a folder
read_wm <- function(path) {
    wm_obj <- list()
    field_names <- c("data_igs", "data_rtrials", "class_labels")
    for (field_name in field_names) {
        wm_obj[[field_name]] <- read_tsv(str_c(path, field_name, ".txt"), progress = FALSE)
    }
    # should be a char vector rather than df
    wm_obj[["class_labels"]] %<>% extract2(1)
    return(wm_obj)
}

# compares clustering quality on differents feature sets and clustering algorithms
watermelon_on_grid <- function(
    data_mat,
    classes,
    feature_sets = list(all_features = rownames(data_mat)),
    # all methods with monotonously increasing distance
    hc_methods = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2"),
    # how much random trials in each case?
    random_trials = 300,
    # how much bootstrap trials for CI calculation?
    bootstrap_trials = 500,
    # display progress bar?
    progress = TRUE,
    # drop a warning if the maxima is non-significant?
    verbose = TRUE,
    # seed for the clusterings
    seed = 451
    ) {

    if (progress) {
        pb <-
            progress_bar$new(
                total = length(hc_methods) * length(feature_sets),
                format = "Exploring grid: [:bar] :percent Time left::eta"
            )
    }
    # holds all the results
    wm_summary <-
        data_frame(
            # all combos of methods and feature sets
            hc_method = rep(hc_methods, each = length(feature_sets)),
            feature_set = rep(names(feature_sets), length(hc_methods)),
            # metric and CI
            wm_area = 0,
            wm_a_upper = 0, 
            wm_a_lower = 0,
            # other things
            random_trials = random_trials,
            bootstrap_trials = bootstrap_trials
        )

    for (hc_method in hc_methods) {
        for(fset_name in names(feature_sets)) {
            fset <- feature_sets[[fset_name]]
            # fast distance; restoring lost dimnames manually
            # plugging it into the clustering
            data_dist <- data_mat[fset, , drop = FALSE] %>% t() %>% parDist() %>% data.matrix()
            dimnames(data_dist) <- list(colnames(data_mat), colnames(data_mat))
            data_dist %<>% as.dist()
            data_hc <- hclust(data_dist, method = hc_method)
            wm <-
                evaluate_hclust(
                    data_hc,
                    classes = classes,
                    random_trials = random_trials,
                    # slightly different for each run
                    seed = seed + which(hc_method == hc_methods) + which(fset_name == names(feature_sets)),
                    progress = FALSE
                )
            wm_boot <- watermelon_bootstrap(wm, step = -1, bootstrap_trials = bootstrap_trials)
            wm_boot_ci <- watermelon_area_ci(wm_boot)
            
            df_mask <-
                wm_summary[["hc_method"]] == hc_method &
                wm_summary[["feature_set"]] == fset_name

            wm_summary[df_mask, ][["wm_area"]] <- watermelon_area(wm)
            wm_summary[df_mask, ][["wm_a_upper"]] <- wm_boot_ci[["upper_bound"]]
            wm_summary[df_mask, ][["wm_a_lower"]] <- wm_boot_ci[["lower_bound"]]
            if(progress) pb$tick()
        }
    }
    # drops a warning if non-significant maxima (CI overlaps with any other point)
    if (verbose) {
        wm_max_i <- which.max(wm_summary[["wm_area"]])
        wm_max_lower <- wm_summary[wm_max_i, ][["wm_a_lower"]]
        wm_rest_upper <- max(wm_summary[-wm_max_i, ][["wm_a_upper"]])
        # if upper bound of any CI criscrosses maxima's lower bound
        if (wm_rest_upper >= wm_max_lower) warning("Non-significant maxima produced!")
    }
    return(wm_summary)
}

plot_grid_result <- function(grid_result, n_best = 10, size = 1) { #, across = c("hc_method", "feature_set")) {
    main_color <- "deepskyblue4"
    if (n_best < nrow(grid_result)) message(str_c("Showing only ", n_best, " best results of total ",  nrow(grid_result)))
    grid_result %<>%
        mutate(
            trial = str_c(hc_method, "\n", feature_set) %>% forcats::fct_reorder(wm_area)
        ) %>%
    # somethimes there's too much combinations
    top_n(n_best, wm_area)

    ggplot(grid_result, aes(y = trial, x = wm_area, xmin = wm_a_lower, xmax = wm_a_upper, color = hc_method)) +
        geom_point(size = size) +
        geom_errorbarh(size = size, height = 0.2) +
        geom_vline(xintercept  = c(0, 1), color = c("red", "forestgreen")) +
        labs(x = "Watermelon area", y = "Method and feature set")
}

heatmap_grid_result <- function(grid_result) {
    grid_result %<>%
        mutate(
            wm_area_best = ifelse(wm_area >= 0.75 * max(wm_area), round(wm_area, 3), "")
        )
    ggplot(grid_result, aes(x = hc_method, y = feature_set, fill = wm_area)) +
        geom_tile() +
        geom_text(aes(label = wm_area_best)) +
        scale_fill_gradientn(
            name = "Watermelon area",
            colours = RColorBrewer::brewer.pal(6, "YlOrRd"),
            limits = 0:1,
            na.value = "grey50"
        ) +
    theme(
        legend.position = "bottom"
    ) +
    guides(fill = guide_colorbar(title.position = "left", barwidth = 13, barheight = 1)) +
    labs(x = "Clustering method", y = "Feature set")
}

prefix = "/mnt/disk1/New_archive/Sh_SEQC_GDC/"

Ps = c('571', '10558', '11154', 'AGL', '75')
NP = length(Ps)

Qs = c("0", "1")
NQ = length(Qs)

methods = array("XYZXYZ", dim = c(NP*NQ))

for ( nq in 1:NQ ) {

   for ( np in 1:NP ) {
   
        methods[np + (nq-1)*NP] = paste(Ps[np], Qs[nq], sep = "_")
        
   }

}

M = length(methods)

bar_row_labels = c("Platform", "Samples")
LAB = c(bar_row_labels, "Ratio")

title_for_method = methods

WM = array(dim=c(25,2*M))

for ( iii in 1:25 ) {

if (iii < 10 ) IFN = paste(prefix, "/reduced/meta_0", toString(iii), ".csv", sep = "")
if (iii >= 10 ) IFN = paste(prefix, "/reduced/meta_", toString(iii), ".csv", sep = "")
meta = read.table(IFN, header = TRUE, sep = ",")

NS = nrow(meta)

D = as.vector(meta[,2])
P = as.vector(meta[,3])
S = as.vector(meta[,1])

UD = unique(D)
UP = unique(P)

NUD = length(UD)
CT_D = 1:NUD
t_rgb = col2rgb("red")
t_hsv = rgb2hsv(t_rgb)     
for ( nud in 1:NUD ) {
     t_hsv[1] = (nud-1)/(NUD+1);
     CT_D[nud] = hsv(t_hsv[1],1,t_hsv[3]);
}

NUP = length(UP)
CT_P = 1:NUP
t_rgb = col2rgb("red")
t_hsv = rgb2hsv(t_rgb)     
for ( nup in 1:NUP ) {
     t_hsv[1] = (nup-1)/(NUP+4);
     CT_P[nup] = hsv(t_hsv[1],0.5,t_hsv[3]);
}

CD = array(dim=c(NS))
CP = array(dim=c(NS))

for ( ns in 1:NS ) {
     d = D[ns]
     index = which (d == UD)
     CD[ns] = CT_D[index]
     p = P[ns]
     index = which (p == UP)
     CP[ns] = CT_P[index]
}

for ( m in 1:M )  {

  print (m)

  if (iii < 10 ) IFN = paste(prefix, "reduced/Sh_", methods[m], "_0", toString(iii), ".csv", sep = "")
  if (iii >= 10 ) IFN = paste(prefix, "reduced/Sh_", methods[m], "_", toString(iii), ".csv", sep = "")
  mas_all = read.table(IFN, header = TRUE, sep = ",")
  
  mas_all = as.matrix(mas_all) 
  mas_all = log10(mas_all+1)
  mas_all = t(mas_all)

  DDD = dist(mas_all)
  par(cex.axis = 1) 
  hc = hclust(DDD)
  hcd = as.dendrogram (hc)
  dend_labels = hcd %>% labels
  nlabels = length(dend_labels)
  bar_color_cancer = array("XYZXYZ", dim = c(nlabels))
  bar_color_platform = array("XYZXYZ", dim = c(nlabels))
  for ( ilabel in 1:nlabels ) {
       sample = dend_labels[ilabel]
       index = which ( sample == S )
       bar_color_cancer[ilabel] = CD[index]
       bar_color_platform[ilabel] = CP[index]
       dend_labels[ilabel]=""
  }     
  if ( iii < 10 ) GFN = paste(prefix, "dend/Sh_", methods[m], "_0", toString(iii),  ".ps", sep = "")
  if ( iii >= 10 ) GFN = paste(prefix, "dend/Sh_", methods[m], "_", toString(iii),  ".ps", sep = "")
  postscript(file=GFN)
  par(cex.axis = 1.5) 
  par(cex.main = 2)
  hcd %>% set("labels", dend_labels) %>% plot()
  par(cex.lab = 1.5)
  par(mar = c(8,4,4,5)+.1)
  title(main = title_for_method[m], sub = "", ylab = "Distance")

  bar_colors = cbind(bar_color_platform, bar_color_cancer)
  colored_bars(
        colors=bar_colors,
        dend=hcd,
        sort_by_labels_order = FALSE,
        rowLabels = bar_row_labels,
        cex.rowLabels = 1.5,
        y_scale = bar_y_scale(hcd),
        y_shift = -bar_y_scale(hcd)
  )
  dev.off()
  
  res1 = evaluate_hclust(tree=hc,classes=D,random_trials = 25)
  res10 = watermelon_area(res1)
  WM[iii,m] = res10
  res2 = evaluate_hclust(tree=hc,classes=P,random_trials = 25)
  res20 = watermelon_area(res2)
  WM[iii,m+M] = res20
  
}  

}

CN = array("XYZXYZ",dim = c(2*M))
CN1 = array("XYZXYZ",dim = c(2*M))

for (kkk in 1:3) {

    for ( m in 1:M ) {
       
        if (kkk < 3) CN[m + (kkk-1)*M] = paste (methods[m], bar_row_labels[kkk], sep = "_")  
        CN1[m + (kkk-1)*M] = paste (methods[m], LAB[kkk], sep = "_")  
    
    }

}

OFN = paste(prefix, "/dend/WM10.csv", sep = "")
write.table(WM,OFN,col.names = CN, row.names = FALSE, sep=",")

WM0 = array(dim=c(25,M))

for ( iii in 1:25 ){
   for ( m in 1:M ) {
     
       WM0[iii,m] = (WM[iii,m]*NUD) / (WM[iii,m+M]*NUP)
       
   }    
}

WM1 = cbind(WM,WM0)

OFN = paste(prefix, "/dend/WM101.csv", sep = "")
write.table(WM1,OFN,col.names = CN1, row.names = FALSE, sep=",")

WM1_ = read.table(OFN, header = TRUE, sep=",")
WM0_ = WM1_[,-(1:6)]

GFN = paste(prefix, "/dend/Ratio10.png", sep = "") 
png(file=GFN)
par(cex.axis=1.5, cex.lab=1.5)
boxplot (WM0)
dev.off()


GFN = paste(prefix, "/dend/legend10.ps", sep = "")
postscript(file=GFN)
plot.new( )
legend("topleft", legend = UD, fill = CT_D, border = CT_D, bty = "n", cex = 1.9)
legend("topright", legend = UP, fill = CT_P, border = CT_P, bty = "n", cex = 3) 
dev.off()
  
