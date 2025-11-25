# -*- coding: utf-8 -*-
# R 3.5
# author: yingxiang.li
# plot: structure

# package
options(warn =-1)

# library(showtext)
# showtext_auto(enable = TRUE)


# const
# RB_CLR <- c("#bfbfbf", "#ff0000", "#99e6ff", "#ff9966", "#404040", "#0099e6", "#ff99e6", "#339933", "#800080",
#      "#00ff00", "#0000ff", "#ff00ff", "#ffe699", "#b24d00", "#00ffff",
#     "#808000", "#ff9999", "#008080", "#99bf26", "#7326e6", "#26bf99", "#808080",
#     "#0d660d", "#ffe6e6",
#     "#993333", "#ff6600", "#33004d")
RB_CLR <- c(
    "#0c5061", "#D95F02", "#7570B3", "#E7298A", "#57725b", "#E6AB02",
    "#A6761D", "#666666", "#3c8fa4", "#FCCDE5", "#80B1D3", "#FB8072", 
    "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
    "#B15928", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5"
)
# RB_CLR <- c("#ff994d", "#0099e6", "#ff99e6", "#339933", "#800080", "#ff004d", "#00ff00", "#0000ff", "#ff00ff", "#ffe699", "#b24d00", "#00ffff", "#808000", "#ff9999", "#008080", "#99bf26", "#7326e6", "#26bf99", "#808080", "#0d660d", "#bfbfbf", "#ff0000", "#99e6ff", "#ff9966", "#404040", "#ffe6e6", "#993333", "#ff6600", "#33004d")


# inp
args <- commandArgs()
inp = args[6]


# main
main <- function(inp){
    struc_list <- c()
    k_list <- c()
    for (name in list.files(dirname(inp))) {
        if (endsWith(name, '.Q')) {
            struc_list <- c(struc_list, name)
            k_list <- c(k_list, as.numeric(tail(strsplit(name, '\\.')[[1]], n=2)[1]))
        }
    }
    out_dir <- paste0(dirname(inp), '/stru.k', min(k_list), '-', max(k_list), '/')
    dir.create(out_dir)
    rslt <- prep_data(k_list, inp, out_dir, 3, 15, T, '', '')

    plot_struc(k_list, out_dir, rslt$pop_list, rslt$ind_num, rslt$add,
        rslt$id_clus, rslt$clr)
    plot_struc(k_list, out_dir, rslt$pop_list, rslt$ind_num, rslt$add,
        rslt$id_clus, rslt$clr, T, bar_hgt=300, pop_sep=50)
    plot_struc(k_list, out_dir, rslt$pop_list, rslt$ind_num, rslt$add,
        rslt$id_clus, rslt$clr, F, T, rslt$struc_q_mean)

    plot_struc_each(k_list, out_dir, rslt)

}


# util
# - custo -
get_avg_pop <- function(inp_file, clus_file, clus_col=1, add_size=T, add_name=T){
    # 获取群体的平均值
    inp <- read.table(inp_file)
    clus <- read.table(clus_file, colClasses='character')[, clus_col]
    pop <- unique(sort(clus))
    k = dim(inp)[2]
    df <- array(dim=c(length(pop), k))
    num <- vector(length=dim(df)[1])
    for (i in 1:length(pop)) {
        for (j in 1:dim(inp)[2]) {
            df[i, j] <- mean(inp[which(clus==pop[i]), j])
            num[i] <- sum(clus==pop[i])
        }
    }
    df <- round(df, 3)
    colnames(df) <- paste(1:k, k, sep='/')
    if (add_size) {
        df <- cbind(num, df)
    }
    if (add_name) {
        df <- cbind(pop, df)
    }
    df
}


prep_data <- function(k_list, inp, out_dir, bar_wd=3, bar_wd_spl=bar_wd*5,
        h_clust=F, cut_line='', sgl_slice=''){
    # 准备数据
#   RB_CLR <- make_pal(k_max+2)
    k_max <- max(k_list)
    k_min <- min(k_list)
    clr <- array(dim=c(k_max, k_max))
    clr[k_min, 1:k_min] <- RB_CLR[1:k_min]
    struc_q <- read.table(paste0(inp, '.', k_min, '.Q'))
    if (k_max > k_min) {
        for (x in (k_min+1):k_max) {
            y <- x
            struc_q_new <- read.table(paste0(inp, '.', x, '.Q'))
            struc_q_cor <- cor(struc_q, struc_q_new)
            z <- vector(length=x-1)
            z[1:length(z)] <- F
            top_cor <- vector(length=x)
            for (i in 1:x) {
                top_cor[i] <- max(struc_q_cor[,i])
            }
            top_cor_ord <- order(top_cor, decreasing=T)
            for (i in top_cor_ord) {
                best_match <- which.max(struc_q_cor[, i])
                if (!z[best_match]) {
                    clr[x, i] <- clr[x-1, best_match]
                    z[best_match] <-T
                } else {
                    clr[x,i] <-  RB_CLR[y]
                    y <- y+1
                }
            }
            RB_CLR <- c(clr[x, 1:x], setdiff(RB_CLR, clr[x, 1:x]))
            struc_q <- struc_q_new
        }
    }
    id_clus <- read.table(paste0(inp, '.fam'), colClasses='character')[, 2:1]
    if (!h_clust) {
        pop_list <- unique(id_clus[, 2])
    } else {
        df <- get_avg_pop(paste0(inp, '.', k_min, '.Q'), paste0(inp, '.fam'))
        struc_q <- df[, 3:(k_min+2)]
        pop_count <- c(as.matrix(df[, 1]))
        class(struc_q) <- "numeric"
        struc_q_mean = t(df)
        if (k_max > k_min) {
            for (k in (k_min+1):k_max) {
                df <- get_avg_pop(paste0(inp, '.', k, '.Q'), paste0(inp, '.fam'))
                struc_q_one <- df[, 3:(k_min+2)]
                class(struc_q_one) <- "numeric"
                struc_q <- cbind(struc_q, struc_q_one)
                struc_q_mean <- rbind(struc_q_mean, t(df[,3:dim(df)[2]]))
            }
        }
        pop_list <- pop_count[hclust(dist(struc_q))$order]
        struc_q_mean_file = paste0(out_dir, '/mean.csv')
        write.table(struc_q_mean, struc_q_mean_file, sep=',', col.names=F,
            quote=F)
    }
    if (cut_line[1] != '') {
        front_line <- intersect(pop_list, cut_line)
        if (sgl_slice == ''){
            pop_list <- c(front_line, setdiff(pop_list, cut_line))
        } else {
            pop_list <- front_line
        }
    }
    ind_num <- dim(id_clus)[1]
    pop_num <- length(pop_list)
    if (sgl_slice != '') {
        pop_num <- length(cut_line)
        ind_num <- 0
        for (i in 1:pop_num) {
            ind_num <- ind_num+sum(id_clus[,2] == cut_line[i])
        }
        ind_num <- 5*ind_num
        print(pop_num)
        print(ind_num)
    }
    add <- 0
    if (cut_line[1] != '') {
        sum_cut <- 0
        for (i in 1:length(cut_line)) {
            sum_cut <- sum(id_clus[,2] == cut_line[i])+sum_cut
            add <- (bar_wd_spl-bar_wd)*sum_cut
        }
    }
    rslt <- list(ind_num=ind_num, add=add, pop_list=pop_list, id_clus=id_clus,
        clr=clr, struc_q_mean=struc_q_mean)
    rslt
}


plot_struc_each <- function(k_list, out_dir, rslt) {
    # 绘图，每个 k 一张图
    k_max <- max(k_list)
    k_min <- min(k_list)
    for (k in k_min:k_max) {
        pref <- paste0(out_dir, 'k', k, '.')
        plot_struc(c(k, k), pref, rslt$pop_list, rslt$ind_num, rslt$add,
            rslt$id_clus, rslt$clr)
        plot_struc(c(k, k), pref, rslt$pop_list, rslt$ind_num, rslt$add,
            rslt$id_clus, rslt$clr, T, bar_hgt=300, pop_sep=50)
        plot_struc(c(k, k), pref, rslt$pop_list, rslt$ind_num, rslt$add,
            rslt$id_clus, rslt$clr, F, T, rslt$struc_q_mean)
    }
}


plot_struc <- function(k_list, graph_pref, pop_list, ind_num, add, id_clus, clr,
        if_mono=F, if_mean_df=F, struc_q_mean='', if_png=F, sgl_slice='',
        cut_line='', pad=200, bar_sep=10, pop_sep=10, bar_wd=3,
        bar_wd_spl=bar_wd*3, bar_hgt=200, k_font=2, lab_font=2) {
    # 绘图
    k_max <- max(k_list)
    k_min <- min(k_list)
    k_num <- c(k_max-k_min+1, 1)[(sgl_slice != '')+1]

    if (graph_pref == '') {
        graph_pref <- paste0(inp, '.admix')
    }
    pop_num <- length(pop_list)
    max_num <- 0

    if (if_mean_df) {
        ind_num <- pop_num
        id_clus <- cbind(struc_q_mean[1, ], struc_q_mean[1, ])
        graph_pref <- paste0(graph_pref, 'mean')
        bar_wd <- bar_wd*10
        bar_hgt <- bar_hgt/2
    } else if (if_mono) {
        max_num <- max(table(id_clus[, 2]))
        graph_pref <- paste0(graph_pref, 'mono')
        bar_wd <- bar_wd
        bar_hgt <- bar_hgt
    } else {
        graph_pref <- paste0(graph_pref, 'norm')
    }

    if (if_png) {
        png(file=paste0(graph_pref, '.png'),
            width=2*pad+(pop_num-1)*pop_sep+ind_num*bar_wd+add,
            height=2*pad+(k_num-1)*bar_sep+k_num*bar_hgt)
    } else {
        pdf(file=paste0(graph_pref, '.pdf'),
            width=(2*pad+(pop_num-1)*pop_sep+ind_num*bar_wd+add)/72*(
                (pop_num <= 10)+1),
            height=(2*pad+(k_num-1)*bar_sep+k_num*bar_hgt)/72)
    }

    plot(rbind(c(0, 0), c(2*pad+(pop_num-1)*pop_sep+c(ind_num, max_num*pop_num)[
        if_mono+1]*bar_wd+add, 2*pad+(k_num-1)*bar_sep+k_num*bar_hgt)),
        type='n', axes=F, xlab=NA, ylab=NA)
    bot <- pad+bar_hgt/2
    left <- pad+1

    for (k in k_min:k_max) {
        if (sgl_slice=='' | sgl_slice==k) {
            text(pad-bar_sep, bot, labels=paste('K=', k, sep=''), adj=1,
                cex=k_font)
            bot <- bot+bar_hgt+bar_sep+1
        }
    }

    for (pop in pop_list) {
        pop_one <- which(id_clus[,2] == pop)
        bar_wd_new <- bar_wd
        if (sum(cut_line == pop) == 1) {
            bar_wd_new <- bar_wd_spl
        }
        if (if_mono) {
            bar_wd_new <- bar_wd_new*max_num/length(pop_one)
        }
        axis_x <- round(left+length(pop_one)*bar_wd_new/2)
        text(axis_x, pad-bar_sep, labels=pop, srt=90, adj=1, cex=lab_font)
        text(axis_x, pad+length(k_min:k_max)*(bar_hgt+bar_sep), labels=pop,
            srt=90, adj=0, cex=lab_font)
        bot <- pad+1
        for (k in k_min:k_max) {
            if (sgl_slice == '' | sgl_slice == k) {
                if (if_mean_df) {
                    mean_k_index = paste(1:k, k, sep='/')
                    struc_q <- apply(t(struc_q_mean[mean_k_index, ]), 2,
                        as.numeric)
                } else {
                    struc_q <- read.table(paste0(inp, '.', k, '.Q'))
                }
                horiz <- vector(length=length(pop_one))
                horiz[1:length(horiz)] <- bot
                for (i in 1:k) {
                    horiz_new <- horiz+c(as.matrix(struc_q[pop_one, i]))*bar_hgt
                    act_one <- which(struc_q[pop_one, i] >= 0.05/bar_hgt)
                    if ((length(act_one) > 0) & (sgl_slice == '' |
                        sgl_slice == k)) {
                            rect(left+((1:length(pop_one)-1)*bar_wd_new)[act_one],
                                horiz[act_one],
                                (left+(1:length(pop_one))*bar_wd_new)[act_one],
                                horiz_new[act_one],
                                col=clr[k, i], border=clr[k, i])
                    }
                    horiz <- horiz_new
                }
                bot <- bot+bar_hgt+bar_sep
            }
        }
        left <- left+pop_sep+bar_wd_new*length(pop_one)
    }
    dev.off()
}


# proc
main(inp)
