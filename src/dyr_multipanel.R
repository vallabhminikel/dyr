# STARTUP #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/dyr')
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(lubridate)
  library(janitor)
  library(openxlsx)
  library(magick)
  library(smoother)
  library(beeswarm)
  library(broom)
  library(drc)
  select = dplyr::select
  summarize = dplyr::summarize
  
})

# Utility helpers (formerly sourced from ../helper.R; inlined for a self-contained repo).
# Note: percent, upper, lower, alpha, ci_alpha, clipcopy, and rbind_files are
# (re)defined later in this script and intentionally override the helper versions.
cast_quietly = function(x) suppressWarnings(as.numeric(x))

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

meansd = function(x, digits=1) {
  paste0(formatC(mean(x, na.rm=T),format='f',digits=digits),'±',formatC(sd(x, na.rm=T),format='f',digits=digits))
}


# OUTPUT STREAMS #### 

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
supplement_directory = tibble(name=character(0), title=character(0))



write_supp_table = function(tbl, title='') {
  
  tbl <- tbl %>%
    mutate(across(where(is.Date), ~ format(., "%Y-%m-%d"))) %>%
    mutate(across(where(lubridate::is.POSIXt), ~ format(as.Date(.), "%Y-%m-%d")))
  
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=3,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/data-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}




# FUNCTIONS ####


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

rank_uniq = function(x) match(x, sort(unique(x)))

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = u95 = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}

lower = l95 = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot


rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}




make_sf = function(survdata) {
  sf = survfit(Surv(dpi, acm) ~ cohort_disp, data=survdata)
  sf$cohort = gsub('cohort_disp=','',names(sf$strata))
  sf$color = meta$color[match(sf$cohort, meta$cohort_disp)]
  return(sf)
}

plot_sf = function(sf, anno_dpi=NULL,xlims=c(0,200),leg=TRUE,xlab='days post-inoculation') {

  xlims <- xlims
  ylims = c(0,1.05)

  plot(sf, col=sf$color, lwd=2, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text=xlab,cex=0.8)
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
  mtext(side=2, line=2.75, text='survival',cex=0.8)
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    
    points(
      x = anno_dpi,
      y = rep(max(ylims), length(anno_dpi)),
      pch = 25,                     
      bg = '#000000',        
      col = '#000000',
      cex = 0.5
    )
    
  }
  
  
  if (leg)  {
    legend <- if (!is.null(sf$strata)) {
      sf$cohort
    } else {
      unique(as.character(na.omit(survdata$cohort)))
    }
    
    
    legend(x='bottomleft', y=1, bty='n', cex=.65, lwd=2, legend=legend,
           col=sf$color,
           text.col=sf$color)
    par(xpd=F)
  }
}


plot_wts = function(weights_indiv, anno_dpi=NULL, xlims,ylims=c(-.22, 0.6),leg=TRUE, euth_line=T) {
  
  rel_weights_smry<- weights_indiv %>%
    group_by(wtdpi, color,cohort_disp) %>%
    summarize(.groups='keep',
              mean = mean(rel),
              sd = sd(rel),
              n = n(),
              l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
              u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
    filter(n > 1) 
  
  xlims = xlims
  
  ylims = ylims
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=seq(-2,15,.2), labels=NA, las=2, tck=-0.05)
  axis(side=2, at=seq(-2,15,.2), labels=percent(seq(-2,15,.2), signed=T), lwd=0, las=2, line=-0.5)
  axis(side=2, at=seq(-2,15,.1), labels=NA, las=2, tck=-0.025)
  abline(h=0)
  if (euth_line) {
    abline(h=-.2, col='red', lty=3)
    text(
      x = (xlims[2])*0.98,     
      y = -0.15,
      labels = "endpoint",
      col = "red",
      adj = c(1, 1),
      cex = 0.8
    ) 
  }
  mtext(side=2, line=2.75, text='weight (% change)',cex=0.8)
  
  for (this_cohort in unique(rel_weights_smry$cohort_disp)) {
    subs = rel_weights_smry %>% filter(cohort_disp==this_cohort) 
    points(subs$wtdpi, subs$mean, type='l', lwd=2, col=subs$color)
    polygon(x=c(subs$wtdpi,rev(subs$wtdpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,ci_alpha),border=NA)
  }
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    
  }
  
  if (leg)  {
    legend <- unique(rel_weights_smry[c("cohort_disp","color")])
    legend <- legend[order(legend$cohort_disp), ]  
    
    
    legend(x='bottomright', y=1, bty='n', cex=.8, lwd=2, legend=legend$cohort_disp,
           col=legend$color,
           text.col=legend$color)
    par(xpd=F)
    
  }
}



plot_weights_max<- function(rel_weights, xlims=c(0, 500),ylims=c(-1,0.5), euth_level=-0.2, ylab='weight to max', legend_pos='bottomright'){
  
  xlims = xlims
  ylims = ylims
  
  
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=seq(-2,15,.2), labels=NA, las=2, tck=-0.05)
  axis(side=2, at=seq(-2,15,.2), labels=percent(seq(-2,15,.2), signed=T), lwd=0, las=2, line=-0.5)
  axis(side=2, at=seq(-2,15,.1), labels=NA, las=2, tck=-0.025)
  abline(h=0)
  abline(h=euth_level, col='red', lty=3)
  text(
    x = xlims[2]-10,     
    y = euth_level+0.05,
    labels = "endpoint",
    col = "red",
    adj = c(1, 1),
    cex = 0.8
  )
  mtext(side=2, line=2.75, text=ylab,cex=0.8)
  for (this_cohort in unique(rel_weights$cohort_disp)) {
    subs = rel_weights %>% filter(cohort_disp==this_cohort) 
    points(subs$wtdpi, subs$mean, type='l', lwd=2, col=subs$color)
    polygon(x=c(subs$wtdpi,rev(subs$wtdpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,ci_alpha),border=NA)
  }
  par(xpd=T)
  
  
  
  legend <- unique(rel_weights[c("cohort_disp","color")])
  legend <- legend[order(legend$cohort_disp), ]  
  
  
  legend(x=legend_pos, y=1, bty='n', cex=.8, lwd=2, legend=legend$cohort_disp,
         col=legend$color,
         text.col=legend$color)
  par(xpd=F)
  
  
  
  
  
  
}









plot_nest <- function(nests, anno_dpi = NULL, xlims,ylims = c(0, 2.1)) {
  
  
  
  if (is.finite(max(nests$nstdpi,na.rm = T)) ) {
    xlims <-c(0, max(nests$nstdpi,na.rm = T)+10)
  } else {
    xlims <-xlims
  }
  ylims = c(0, 2.1)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
  
  axis(side=1, at=seq(0,600,100), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,600,50), lwd=0, labels=seq(0,600,50), line=-0.5)
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation')
  axis(side=2, at=seq(0,2,0.5), lwd=0, labels=seq(0,2,0.5), line=-0.5, las=2)
  axis(side=2, at=seq(0,2,0.5), labels=NA, tck=-0.025)
  mtext(side=2, line=2.75, text='nest (mean)')
  
  # thin line for each cage
  for (current_cage in unique(nests$cage)) {
    subs = nests %>% filter(cage==current_cage) 
    points(subs$nstdpi, subs$mean, type='l', col=subs$color[1],lwd=0.25)
  }
  
  # trend line for each group
  for (this_cohort in unique(nests$cohort_disp)) {
    subs = nests %>% 
      filter(cohort_disp==this_cohort) %>% 
      group_by(cohort_disp,nstdpi,color) %>% 
      summarize(.groups='keep',
                mean=mean(mean)) %>%
      ungroup()
    gaussian_smooth <- smth(subs$mean, method = "gaussian", window = 5, tails = TRUE)
    points(subs$nstdpi, gaussian_smooth, type = 'l', col = subs$color[1], lwd = 2)
    
    
  }
  
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    points(x=c(anno_dpi), y=rep(max(ylims),length(anno_dpi)), pch=25, bg='#000000',cex=0.15)
  }
  
  
  legend <- unique(nests[c("cohort_disp","color")])
  legend <- legend[order(legend$cohort_disp), ]  
  
  
  
  
  legend(x='bottomleft', y=1, bty='n', cex=.8, lwd=2, legend=legend$cohort_disp,
         col=legend$color,
         text.col=legend$color)
  
  par(xpd=F)
  
  
}




intmode = function(x) {
  as.integer(table(x)[which.max(table(x))])
}


smooth_p = function(p) {
  if (length(p)<2) {
    return (p)
  }
  while (any(cummin(p) != p)) { # while p is non-monotonic
    for (i in 2:length(p)) { # starting from the second dilution and going up,
      if (p[i] > p[i-1]) { # if the current one has higher p than the last one...
        # select all dilutions below the current one that have l95 p...
        indices_to_average = c(which(1:length(p) < i & p < p[i]), i)
        # and average them all together with the current value.
        p[indices_to_average] = mean(p[indices_to_average])
      }
    }
  }
  return (p)
}



spearman_karber = function(x, n, d) {
  # first sort everything from strongest to dilutest
  sorted = tibble(x=x, n=n, d=d) %>%
    arrange(desc(d))
  x = sorted$x
  n = sorted$n
  d = sorted$d
  # to formally define LLQ and ULQ, assume all samples would be positive at one dilution higher
  # and all negative at one dilution l95
  if (length(x)==1) {
    dd = 1
  } else {
    dd = mean(abs(diff(d))) # find dd, the difference between dilutions
    if (any(abs(diff(d)) - dd > .05)) { # .05 provides some numeric tolerance, esp in CSF dilution where rounding error is substantial
      stop("x are not uniformly distributed\n")
    }
  }
  n_mode = intmode(n) # the most common number of replicates
  x = c(n_mode, x, 0) # to define llq & ulq, assume it'd be 100% pos at 1 higher dilution and 0% pos at 1 l95 dilution
  n = c(n_mode, n, n_mode) # assume the higher and l95 dilution would've been the modal N
  d = c(d[1]+dd, d, d[length(d)]-dd) # assume the higher and l95 dilution would have followed the same dilution factor
  llq = mean(d[1:2]) # define the LLQ as what you'd get if you had 100% pos at 1 higher dilution but 0% at highest dilution tested
  ulq = mean(d[(length(d)-1):length(d)]) # define the ULQ as what you'd get if you had 0% pos at 1 l95 dilution but 100% at lowest dilution tested
  # if a non-zero positive is observed at a dilution higher than where zeroes are observed, call this noise and remove it
  for (j in 2:length(x)) {
    if (x[j-1]==0 & x[j] > 0) {
      x[j] = 0
    }
  }
  p = smooth_p(x/n) # smooth the proportion positive to be non-increasing
  nd = length(x) # number of dilutions
  # find i, the dilutest dilution at which 100% positive replicates
  if (all(p < 1)) {
    i = 1
  }  
  else  {
    i = max(which(p==1))
  }
  sd50 = (d[i]+dd) - dd*(1/2 + sum(p[i:nd])) # calculate log10(sd50)
  se = sqrt(dd^2 * sum(p*(1-p)/(n-1))) # calculate standard error
  # invert estimate and limits of quantification, in order to be in titer space rather than dilution space:
  return_value = cbind(estimate=-sd50, se=se, llq=-llq, ulq=-ulq)
  return (return_value)
}







# DATA READ IN  ####
nest= read_tsv('data/analytic/nests.tsv', col_types=cols())
weight= read_tsv('data/analytic/weights.tsv', col_types=cols())
blind = read_tsv('data/analytic/blind.tsv', col_types=cols())
master=read_tsv( 'data/analytic/master.tsv', col_types=cols())
surgery =read_tsv( 'data/analytic/surgeries.tsv', col_types=cols())
ellas = read_tsv('data/analytic/dyr_ellas_broad.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(mean_conc = suppressWarnings(as.numeric(gsub('<','',mean_conc))))
elisas = read_tsv('data/analytic/dyr_elisas.tsv', col_types=cols())
gt_meta = tibble(genotype=c('WT','HE','KO'),
                 color=c('#0001CD','#AB23AB','#AD2016'))

rat_it_icv = read_tsv('data/analytic/ionis_rat_it_vs_icv.tsv', col_types=cols())

## PREP ####

meta = read_tsv('data/analytic/meta.tsv', col_types=cols()) %>%
  distinct(expt,cohort,cohort_disp,species,genotype,inoculum,treatment,dose_ug,treatment_start_dpi,color)

cohort_lookup <- read_tsv('data/analytic/cohort_lookup.tsv', col_types=cols())%>% clean_names()%>%
  rename(
    cohort = `cohort_name`,
    display = `display_cohort_name`
  )

### methods ####


master %>%
  inner_join(blind %>% select(-cage), by='animal') %>%
  filter(!(inoculum %in% c('none','uninoculated'))) %>%
  mutate(age_days_at_inoc = as.integer(inoculation_date - dob)) %>%
  group_by(species) %>%
  summarize(.groups='keep',
            n=n(),
            min_age_weeks = min(age_days_at_inoc, na.rm=T)/7,
            max_age_weeks = max(age_days_at_inoc, na.rm=T)/7) %>%
  ungroup() -> age_range_at_inoc

write_supp_table(age_range_at_inoc, 'Age ranges of animals at inoculation by species.')


### surv ####



survdata<- master %>% 
  inner_join(select(blind, species,expt, animal, cohort,inoculum,genotype,treatment ,treatment_start_dpi), by=c('animal','expt')) %>%
  mutate(death_date = as.Date(death_date),
         prion_endpoint=as.logical(prion_endpoint)) %>%
  select(expt,animal, species,cohort, genotype,txdpi=treatment_start_dpi, treatment, dob, inoculum,inoculation_date,death_date, prion_endpoint, acm,endpoint_comments,acm_comments) %>%
  filter(!str_detect(cohort, "-died|excluded"))%>%
  inner_join(meta %>% distinct(expt,color,cohort,cohort_disp), by=c('expt','cohort')) %>%
  mutate(dpi = as.integer(death_date - inoculation_date))

survdata$txdpi = as.character(unlist(survdata$txdpi))
survdata$acm=as.logical(survdata$acm)

###weight ####


weight_plot <- weight %>%
  mutate(weight_date = as.Date(weight_date)) %>%
  inner_join(master %>% 
               select(expt, animal, inoculation_date, prion_endpoint, acm),by=c('animal','expt') )%>%
  mutate(inoculation_date = as.Date(inoculation_date)) %>%
  inner_join(blind %>% select(expt,cohort,animal,txdpi=treatment_start_dpi), by=c('animal','expt'))%>%
  inner_join(meta %>% distinct(expt,color,cohort,cohort_disp), by=c('expt','cohort'))%>%
  filter(!str_detect(cohort, "-died|excluded"))%>%
  mutate(wtdpi = as.integer(weight_date - inoculation_date)) %>%
  group_by(expt, animal,cohort_disp) %>%
  mutate(
    baseline = first(weight),
    rel = (weight - baseline) / baseline,
    peak     = max(weight, na.rm = TRUE)   ,
  ) %>%
  ungroup()%>%
  mutate(offpeak = (weight - peak) / peak) 



### nest ####

nest_plot <- nest %>%
  mutate(nest_date = as.Date(nest_date),
         cage=unlist(cage)) %>%
  left_join(master %>%
              select(expt,cage, inoculation_date,animal),by=c('cage','expt'), relationship='many-to-many')%>%
  mutate(inoculation_date = as.Date(inoculation_date)) %>%
  left_join(blind %>% 
              distinct(animal,cohort,treatment,cage,txdpi=treatment_start_dpi,genotype), by=c('cage','animal'))%>%
  inner_join(meta %>% distinct(expt,color,cohort,cohort_disp), by=c('expt','cohort'))%>%
  group_by(cage,expt) %>%
  filter(n_distinct(genotype, na.rm = TRUE) == 1) %>%  
  ungroup()%>%
  filter(!is.na(inoculation_date)) %>%
  filter(!str_detect(cohort, "-died|excluded"))%>%
  mutate(nstdpi = as.integer(nest_date - inoculation_date))%>%
  rename (tx=treatment, comb=combined_score) %>%
  mutate(mean=comb)%>%
  select(-comb)


stat_registry <- new.env()
stat_registry$surv <- list()
stat_registry$peak <- list()

describe_groups <- function(dat) {
  has_mri <- exists('meta_mri') && exists('blind_mri')
  mri_expts <- if (has_mri) meta_mri$expt else character(0)
  info <- bind_rows(
    meta %>%
      transmute(expt, cohort, cohort_disp = as.character(cohort_disp), species,
                genotype = recode(genotype, HE = 'HET', KO = 'HOM'),
                prion_strain = as.character(inoculum), treatment = as.character(treatment),
                dose_ug = as.numeric(dose_ug),
                treatment_start_dpi = as.character(treatment_start_dpi)),
    if (has_mri) meta_mri %>%
      left_join(blind_mri %>% distinct(expt, species), by = 'expt') %>%
      transmute(expt, cohort, cohort_disp = as.character(cohort_disp), species,
                genotype = recode(genotype, 'TgVole' = 'Tg(Vole)', 'Prp null' = 'PrP KO', 'Prnp Null' = 'PrP KO'),
                prion_strain = coalesce(as.character(inoculum), 'none'), treatment = as.character(treatment),
                dose_ug = suppressWarnings(as.numeric(dose_ug)),
                treatment_start_dpi = NA_character_)
  )
  present <- if ('cohort' %in% names(dat)) {
    dat %>% distinct(expt, cohort) %>% inner_join(info, by = c('expt', 'cohort'))
  } else {
    dat %>% distinct(expt, cohort_disp) %>% mutate(cohort_disp = as.character(cohort_disp)) %>%
      inner_join(info, by = c('expt', 'cohort_disp'))
  }
  collapse <- function(x) paste(sort(unique(na.omit(as.character(x)))), collapse = '; ')
  present %>%
    group_by(group = cohort_disp) %>%
    summarize(.groups = 'drop',
              experiment = collapse(expt),
              species = collapse(species),
              genotype = collapse(genotype),
              prion_strain = collapse(prion_strain),
              treatment = collapse(treatment),
              dose_ug = collapse(dose_ug),
              treatment_start_dpi = collapse(treatment_start_dpi),
              time_axis = ifelse(all(expt %in% setdiff(mri_expts, 'HU2')), 'age (days)', 'dpi'))
}

surv_stats <- function(dat, ref_group) {
  dat <- dat %>% mutate(group = as.character(cohort_disp))
  sf <- survfit(Surv(dpi, acm) ~ group, data = dat)
  med <- tibble(group = sub('group=', '', names(sf$strata)),
                median_survival = unname(summary(sf)$table[, 'median']))
  pvals <- map_dfr(setdiff(unique(dat$group), ref_group), function(g) {
    dat2 <- dat %>% filter(group %in% c(ref_group, g)) %>%
      mutate(tmp = factor(group, levels = c(ref_group, g)))
    s <- suppressWarnings(survdiff(Surv(dpi, acm) ~ tmp, data = dat2))
    tibble(group = g, p_vs_control = signif(1 - pchisq(s$chisq, df = 1), 3))
  })
  dat %>%
    group_by(group) %>%
    summarize(.groups = 'drop', n = n(), events = sum(acm == 1)) %>%
    left_join(med, by = 'group') %>%
    left_join(pvals, by = 'group') %>%
    mutate(control_group = ref_group)
}

surv_table <- function(dat, ref_group) {

  sf <- survfit(Surv(dpi, acm) ~ cohort_disp, data = dat)
  
  # survfit medians
  km_med <- as.data.frame(summary(sf)$table) %>%
    tibble::rownames_to_column("cohort_disp") %>%
    mutate(
      cohort_disp = sub("cohort_disp=", "", cohort_disp),
      
    ) %>%
    select(cohort_disp, median)
  
  
  # N, events, mean±sd[acm] 
  basic <- dat %>%
    group_by(cohort_disp) %>%
    summarise(
      N = n(),
      experiments = paste(unique(expt), collapse = ", "),
      events = sum(acm == 1),
      mean_acm = mean(dpi[acm == 1], na.rm = TRUE),
      sd_acm = sd(dpi[acm == 1], na.rm = TRUE),
      .groups = "drop"
    )%>% ungroup()
  
  # reference group median
  ref_median <- km_med$median[km_med$cohort_disp == ref_group]
  
  # percent delta from reference
  km_med <- km_med %>%
    mutate(delta = (median - ref_median) / ref_median)
  
  all_groups <- sort(unique(surv_ex$cohort_disp))
  compare_groups <- setdiff(all_groups, ref_group)
  
  pairwise <- map_dfr(compare_groups, function(g) {
    
    dat2 <- surv_ex %>% 
      filter(cohort_disp %in% c(ref_group, g)) %>%
      mutate(tmp = factor(cohort_disp, levels = c(ref_group, g)))
    
    s <- suppressWarnings(survdiff(Surv(dpi, acm) ~ tmp, data = dat2))
    p <- 1 - pchisq(s$chisq, df = 1)
    
    tibble(cohort_disp = g, p_vs_ref = p)
  }) %>%
    bind_rows(tibble(cohort_disp = ref_group, p_vs_ref = 1))
  
  
  final <- basic %>%
    left_join(km_med, by = "cohort_disp") %>%
    left_join(pairwise, by = "cohort_disp") %>%
    mutate(
      median = round(median, 1),
      delta = round(delta, 2),
      mean_acm = round(mean_acm, 1),
      sd_acm = round(sd_acm, 2),
      p_vs_ref = signif(p_vs_ref, 3)
    )%>%   mutate( median = ifelse(
      is.na(median),
      paste0(">", max(dat$dpi, na.rm = TRUE)), # if no median reched then use last dpi
      median
    ))%>%
    mutate(across(everything(), ~ ifelse(is.na(.x) | is.nan(.x), "", .x)))

  stat_registry$surv[[length(stat_registry$surv) + 1]] <- surv_stats(dat, ref_group) %>%
    left_join(describe_groups(dat), by = 'group')
  
  return(final)
}


make_weight_summary <- function(data, ref_group, stat_ref_group = ref_group) {
  
  #  Peak 
  peak_by_animal <- data %>%
    group_by(cohort_disp, animal) %>%
    summarise(
      peak_weight = weight[which.max(weight)],
      peak_rel    = rel[which.max(weight)],       # % change from baseline
      peak_dpi    = wtdpi[which.max(weight)],     # dpi of peak
      .groups = "drop"
    )%>% ungroup()
  
  stat_peaks <- data %>%
    group_by(group = as.character(cohort_disp), animal) %>%
    summarise(peak_dpi = wtdpi[which.max(weight)], .groups = 'drop')
  ref_dpi_vals <- stat_peaks$peak_dpi[stat_peaks$group == stat_ref_group]
  stat_registry$peak[[length(stat_registry$peak) + 1]] <- stat_peaks %>%
    group_by(group) %>%
    summarize(.groups = 'drop',
              n = n(),
              mean_peak_dpi = mean(peak_dpi),
              l95 = mean(peak_dpi) - 1.96 * sd(peak_dpi) / sqrt(n()),
              u95 = mean(peak_dpi) + 1.96 * sd(peak_dpi) / sqrt(n()),
              p_vs_control = ifelse(group[1] == stat_ref_group, NA_real_,
                                    tryCatch(signif(t.test(peak_dpi, ref_dpi_vals)$p.value, 3),
                                             error = function(e) NA_real_))) %>%
    mutate(control_group = stat_ref_group) %>%
    left_join(describe_groups(data), by = 'group')
  
  # Cohort-level peaks
  cohort_peaks <- peak_by_animal %>%
    group_by(cohort_disp) %>%
    summarise(
      n=n(),
      mean_peak_weight = mean(peak_weight, na.rm = TRUE),
      mean_peak_pct    = mean(peak_rel * 100, na.rm = TRUE),
      mean_peak_dpi    = mean(peak_dpi, na.rm = TRUE),
      .groups = "drop"
    )%>% ungroup()
  
  #  Add experiment list per cohort 
  cohort_expts <- data %>%
    group_by(cohort_disp) %>%
    summarise(
      experiments = paste(sort(unique(expt)), collapse = ", "),
      .groups = "drop"
    )%>% ungroup()
  
  
  cohort_peaks <- cohort_peaks %>%
    left_join(cohort_expts, by = "cohort_disp") %>%
    relocate(experiments, .after = cohort_disp)
  
  
  ref_dpi <- cohort_peaks$mean_peak_dpi[
    cohort_peaks$cohort_disp == ref_group
  ]
  
  cohort_peaks <- cohort_peaks %>%
    mutate(
      delta_dpi = (mean_peak_dpi - ref_dpi) / ref_dpi
    )
  
  
  cohort_peaks %>%
    mutate(
      mean_peak_weight = round(mean_peak_weight, 2),
      mean_peak_pct    = round(mean_peak_pct, 2),
      mean_peak_dpi    = round(mean_peak_dpi, 2),
      delta_dpi        = round(delta_dpi, 2)
    ) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x) | is.nan(.x), "", .x)))
}



## FIGURE ####







# RAT ####
## ELISA ####





smry_132 <- read_tsv('data/analytic/elisa_rat_132.tsv', col_types=cols(sample=col_character()))
smry_132$sample=unlist(smry_132$sample)

smry_168 <- read_tsv('data/analytic/elisa_rat_168.tsv', col_types=cols(sample=col_character()))
smry_168$sample=unlist(smry_168$sample)



## save blind for future

blind_hamster<- read_tsv('data/analytic/blind_hamster_132_saves.tsv', col_types=cols())%>%
  mutate(sample=as.character(sample))

blind_rat<- read_tsv('data/analytic/blind_rat_168_save.tsv', col_types=cols())

ana_132<- smry_132%>%
  inner_join(blind_hamster,by='sample')


ana_168<- smry_168%>%
  inner_join(blind_rat,by='sample')%>%
  mutate(genotype=case_when(
    str_detect(genotype, "HO") ~ "HOM",
    str_detect(genotype, "HE") ~ "HET",
    TRUE ~ genotype
  ))






plot_csELISA <- function(df, species,cex_a=0.8, ref_lines=FALSE, geno_labels=c("WT","HET","HOM")) {
  
  
  
  df <- df %>%
    filter(genotype %in% c("WT", "HET", "HOM"))%>%
    mutate(col = case_when(
      genotype == "WT"  ~ "#0001CD",
      genotype == "HET" ~ "#501378",
      genotype == "HOM" ~ "#A02422"))
  
  wt_mean <- mean(df$ngml_av[df$genotype == "WT"], na.rm = TRUE)
  df <- df %>%
    mutate(norm_conc = ngml_av / wt_mean )
  
  
  LLQ <- round(0.05 * 200 / wt_mean, 5) # 0.05 ng/mL standard, 1:200 dilution, norm'ed to WT mean as 100%
  
  n_by_group <- df %>% count(genotype, name = "n")
  x_label <- setNames(
    paste0(geno_labels,
           "\n(n=", n_by_group$n[match(c("WT","HET","HOM"), n_by_group$genotype)], ")"),
    c("WT","HET","HOM")
  )
  
  df$genotype <- factor(df$genotype, levels = c("WT","HET","HOM"))
  
  
  plot(
    NA, NA,
    xlim = c(0.5, 3.5),
    ylim = c(0, max(df$norm_conc) * 1.2),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  
  
  for (i in 1:3) axis(1, at = i, labels = x_label[c("WT","HET","HOM")][i], line = 0, col = NA,  tick = FALSE,cex.axis=cex_a )
  axis(1, at=0:4, labels=NA, tck=FALSE)
  
  axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2)
  par(xpd=T)
  stat <- df %>%
    group_by(genotype,col) %>%
    summarise(
      mean = mean(norm_conc, na.rm = TRUE),
      se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(genotype, levels = c("WT","HET","HOM")))
    )
  
  

  suppressWarnings(with(stat, {
    arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col = col)
    segments(x - 0.25, mean, x + 0.25, mean, lwd = 3,col = col)
  }))


  beeswarm(df$norm_conc ~ df$genotype,
           pch = 21, cex = 0.8,
           pwcol =df$col,
           add = TRUE)
  
  
  text(
    x = par("usr")[2]- 0.5, 
    
    y = LLQ+0.1,
    labels = "LLQ",
    
    col = "red", pos = 4, xpd = TRUE
  ) 
  par(xpd=F)
  
  
  if (ref_lines) {
    abline(h = 1,   lty = 2, col = "gray50")
    abline(h = 0.5, lty = 2, col = "gray80")
  }
  
  abline(h = LLQ, lty = 2, col = "red")
  
  
  
  
  
  return(stat)
  
}


rat_elisa <- function(){ 
  par(mar = c(4,5,2,0))
  
  plot_csELISA(ana_168, "Rat",cex_a=1, ref_lines=TRUE, geno_labels=c("WT","Het KO","Hom KO"))
  
  mtext(side=2, line=3.2, text='residual PrP (%)',cex=0.8)
  
}

stat_rat <- {
  df <- ana_168 %>%
    filter(genotype %in% c("WT", "HET", "HOM")) %>%
    mutate(col = case_when(genotype == "WT"  ~ "#0001CD",
                           genotype == "HET" ~ "#501378",
                           genotype == "HOM" ~ "#A02422"))
  wt_mean <- mean(df$ngml_av[df$genotype == "WT"], na.rm = TRUE)
  df %>%
    mutate(norm_conc = ngml_av / wt_mean) %>%
    group_by(genotype, col) %>%
    summarise(mean = mean(norm_conc, na.rm = TRUE),
              se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    ungroup() %>%
    mutate(l95 = mean - 1.96 * se,
           u95 = mean + 1.96 * se,
           x = as.numeric(factor(genotype, levels = c("WT","HET","HOM"))))
}
write_supp_table(stat_rat %>% select(-x,-col), "Rat ELISA summarized")

## RR  ####

### survival ####

surv_ex <- survdata %>%
  filter(expt %in% c("RR"))


write_supp_table(surv_ex %>% select(-color,-cohort), 'RR survival individual')

sf = make_sf(surv_ex)

surv_smry<- surv_table(surv_ex, ref_group = "NBH")%>%
  arrange(match(cohort_disp, c('NBH','RRML')))

write_supp_table(surv_smry, 'RR survival summarized')



RR_survival <- function(){ 
  par(mar=c(4,5,2,4))
  
  par(xpd=F)
  
  
  surv_ex <- survdata %>%
    filter(expt %in% c("RR"))
  
  sf = make_sf(surv_ex)
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
}







### weight ####

dat_ex<- weight_plot %>%
  filter(expt == 'RR') %>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  

write_supp_table(dat_ex%>%select(-color,-cohort), 'RR weights individual')

rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup() %>%
  filter(n > 1) %>%
  select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'RR weights summarized')


weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt == 'RR'), ref_group = "RRML", stat_ref_group = "NBH")%>%
  arrange(match(cohort_disp, c('NBH','RRML')))

write_supp_table(weight_summary, 'RR weights peak')

RR_weight <- function(){ 
  
  par(mar=c(4,3,2,2))
  
  par(xpd=F)
  
  
  dat_ex<- weight_plot %>%
    filter(expt == 'RR') %>%
    select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  
  
  
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(120, max(dat_ex$wtdpi)),ylims = c(-0.2,0.25),leg=FALSE,euth_line=F)
  
  
  
  
  
  
}






### nfl ####

nfl_rr = ellas %>%
  filter(experiment=='RR',
         is.na(exclude))%>%
  select(animal=animal_id,species,expt=experiment,tissue,sample_dpi,sample_name,mean_rfu,rfu_percent_cv,mean_conc,conc_percent_cv)


plot_nfl_rr<- nfl_rr %>%
  filter(tissue == "plasma") %>%
  left_join(blind %>% select(expt,cohort,animal,genotype), by=c('expt','animal'))%>%
  left_join(meta %>% distinct(expt,cohort,color,cohort_disp), by=c('expt','cohort'))%>%
  group_by(expt,cohort, sample_dpi,color,cohort_disp) %>%
  summarise(
    .groups   = "keep",
    mean = mean(mean_conc, na.rm = TRUE),
    n = n(),
    l95 = lower(mean_conc),
    u95 = upper(mean_conc)
    
  )%>% ungroup()

write_supp_table(plot_nfl_rr%>%select(-color), 'RR NfL summarized')

RR_nfl <- function(){ 
  
  
  par(mar=c(3,5,2,4))
  par(xpd=F)
  
  
  xlims <-c(-10,245)
  ylims =c(0,260)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
  
  axis(side=1, at=seq(-10,240,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,240,50), lwd=0, labels=seq(0,240,50), line=-0.5)
  axis(side=1, at=seq(0,240,50), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=seq(0,250,50), lwd=0, labels=seq(0,250,50), line=-0.5, las=2)
  axis(side=2, at=seq(0,250,10), labels=NA, tck=-0.025)
  mtext(side=2, line=2.5, text='plasma NfL (pg/ml)',cex=0.8)
  
  
  nfl_rr %>%
    filter(tissue == "plasma") %>%
    left_join(blind %>% select(expt,cohort,animal,genotype), by=c('expt','animal'))%>%
    left_join(meta %>% distinct(expt,cohort,color,cohort_disp), by=c('expt','cohort'))->indiv
  par(xpd=T)
  
  points(indiv$sample_dpi,indiv$mean_conc,col=indiv$color)
  par(xpd=F)
  
  for (co in unique(plot_nfl_rr$cohort)) {
    subs <- plot_nfl_rr %>%
      filter(cohort==co)
    
    polygon(x=c(subs$sample_dpi,rev(subs$sample_dpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,0.5),border=NA)
    lines(subs$sample_dpi, subs$mean, lwd = 2.5, col = subs$color)
    points(subs$sample_dpi, subs$mean, pch = 16, col = subs$color)
    
  }
  
  legend(x='topleft', y=1, bty='n', cex=.8, lwd=2, legend=unique(plot_nfl_rr$cohort_disp),
         col=unique(plot_nfl_rr$color),
         text.col=unique(plot_nfl_rr$color))

  ks_by_dpi_rr <- indiv %>%
    filter(cohort_disp %in% c("RRML", "NBH")) %>%
    group_by(sample_dpi) %>%
    summarize(
      n_RRML  = sum(cohort_disp == "RRML"),
      n_NBH   = sum(cohort_disp == "NBH"),
      p_value = if (n_RRML >= 2 & n_NBH >= 2) {
        tryCatch(
          ks.test(
            mean_conc[cohort_disp == "RRML"],
            mean_conc[cohort_disp == "NBH"]
          )$p.value,
          error = function(e) NA
        )
      } else {
        NA
      },
      .groups = "drop"
    ) %>% ungroup() %>%
    mutate(
      sig = case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )

  mtext(at = ks_by_dpi_rr$sample_dpi, line=0.25, text = ks_by_dpi_rr$sig, cex = 0.7, font = 1)


}




### rt-quic ####

# Eric - hamster CSF
csf_unit_factor = 0.001 # the "dilution" col for CSF RT-QuIC reflects the number of microliters of pure CSF equivalent per well. thus, a factor of 0.001 gets you to units of seeds per mL
bh_unit_factor = 0.002 # the "dilution" col for brain homogenate RT-QuIC reflects the final wt/vol homogentae, where the seed is 2 µL (i.e. equivalent 2 mg tissue if 100% wt/vol) so a factor of 0.002 gets you to units of seeds per gram

dyr_rtquics = read_tsv('data/analytic/dyr_rtquics.tsv', col_types=cols()) %>%
  filter((is.na(exclude) | exclude != 'y'))
master = read_tsv('data/analytic/master.tsv', col_types=cols())
blind = read_tsv('data/analytic/blind.tsv', col_types=cols())


rt_quic_rr = dyr_rtquics %>%
  filter(experiment=='RR',
         is.na(exclude),
         tissue=='brain')%>%
  clean_names()%>% 
  mutate(animal_id=unlist(animal_id),
         sample=unlist(sample))%>%
  rename(animal=animal_id,expt=experiment) %>%
  arrange(desc(rta), animal, dilution) %>%
  group_by(animal, dilution) %>% slice(1) %>% ungroup() %>% # for samples re-run, remove the earlier version
  arrange(rta, animal, desc(dilution))


plot_rt_quic_rr <-  rt_quic_rr %>%
  left_join(blind %>% select(animal, cohort, genotype, treatment, inoculum), by=c('animal'))%>%
  left_join(meta %>% select(cohort, cohort_disp,color), by=c('cohort'))%>%
  left_join(master %>%
              select(-cage) %>%
              mutate(dpi = as.integer(death_date - inoculation_date)), by=c('animal','expt'))%>%
  select(expt, species, animal, inoculum, genotype, color, dpi, treatment, dilution, positive_replicates, total_replicates) 


plot_rt_quic_rr %>%
  group_by(expt, species, animal, inoculum, genotype, dpi, treatment, color) %>%
  summarize(.groups='keep',
            sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution)))%>% ungroup() %>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2]) %>%
  select(-sk) %>% ungroup() -> titers_rr

write_supp_table(titers_rr%>%select(-color), 'RR RT-quic summarized')


titers_rr$treatment <- factor(titers_rr$inoculum, levels = c("NBH", "RRML"))


RR_rtqic <- function(){ 
  par(xpd=F)
  
  par(mar=c(3,3,2,2))
  
  plot(
    NA, NA,
    xlim = c(0.5, 2.5),
    ylim = c(4,12),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  
  axis(1, at = 1:2, labels = c("NBH", "RRML"), line = -1, col = NA,  tick = FALSE )
  axis(1, at=0:3, labels=NA, tck=FALSE)
  
  axis(2, las = 2)
  axis(1,labels = FALSE,tick = FALSE)
  mtext(side=2, line=2, text='brain prion titer\nlog10 seeds/g',cex=0.8)
  
  stat <- titers_rr %>%
    group_by(inoculum,color) %>%
    summarise(
      mean = mean(titer, na.rm = TRUE),
      se = sd(titer, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(inoculum, levels = c("NBH", "RRML")))
    )
  
  
  
  LLQ <- titers_rr %>% 
    filter(titer > 0) %>% 
    summarise(llq = min(titer, na.rm = TRUE)) %>% 
    pull(llq)
  
  
  abline(h = LLQ, lty = 2, col = "red")
  
  
  
  mtext(side=4, line=0.1, at=5.69, text="LLQ", col="red", las=2, cex=0.7)
  
  
  with(stat, {
    suppressWarnings(arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col = color))
    segments(x - 0.25, mean, x + 0.25, mean, lwd = 3,col = color)
  })
  
  
  beeswarm(titers_rr$titer ~ titers_rr$inoculum,
           pch = 21, cex = 1.2,
           pwcol =titers_rr$color,
           add = TRUE)
  
  
  
  
  
}

write_supp_table( titers_rr %>%
                    group_by(inoculum,color) %>%
                    summarise(
                      .groups = "keep",
                      mean = mean(titer),
                      l95 = lower(titer),
                      u95 = upper(titer)
                    ) %>% ungroup()%>%
                    mutate(  
                      x = as.numeric(factor(inoculum, levels = c("NBH", "RRML")))
                    )%>%select(-color,-x), 'RR RT-quic summary')



### prp ELISA -- data pending ####













## RD ####

### CSF PrP ####



elisa_rd = elisas %>%
  filter(is.na(exclude),
         experiment=='RD')%>%
  mutate(animal_id=as.character(unlist(animal_id)),
         sample=unlist(sample),
         ngml_av=unlist(ngml_av),
         se_mean=unlist(se_mean))%>%
  rename(animal=animal_id,expt=experiment)

plot_elisa_rd <- elisa_rd %>%
  left_join(blind %>% select(animal, cohort, genotype, treatment, inoculum), by=c('animal'))%>%
  mutate(col = case_when(
    genotype == "WT" & treatment=='ASO 6' ~ "#0001CD",
    genotype == "HET" ~ "#501378",
    genotype == "HOM" ~ "#A02422",
    treatment=='PBS' ~ "#A9A9A9"))







aso_mean <- mean(plot_elisa_rd$ngml_av[plot_elisa_rd$treatment == "PBS"], na.rm = TRUE)

plot_elisa_rd <- plot_elisa_rd %>%
  mutate(norm_conc = ngml_av / aso_mean )


write_supp_table(plot_elisa_rd %>% select(-exclude,-col,-cohort), 'RD CSF PrP individual')


LLQ <- 10 / aso_mean * 100

plot_elisa_rd %>%
  group_by(treatment,col) %>%
  summarise(
    mean = mean(norm_conc, na.rm = TRUE),
    se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>% ungroup()%>%
  mutate(
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
    x = as.numeric(factor(treatment, levels = c("PBS", "ASO 6")))
  ) ->stat_save

write_supp_table(stat_save%>% select(-col,-x)%>%
                   arrange(match(treatment, c('PBS','ASO 6'))), 'RD CSF PrP summarized')





RD_CSF<- function(){ 
  
  
  n_by_group <- plot_elisa_rd %>% count(treatment, name = "n")
  x_label <- setNames(
    paste0(levels(factor(plot_elisa_rd$treatment, levels = c("PBS", "ASO 6"))),
           "\n(n=", n_by_group$n[match(c("PBS", "ASO 6"), n_by_group$treatment)], ")"),
    c("PBS", "ASO 6")
  )
  
  plot_elisa_rd$treatment <- factor(plot_elisa_rd$treatment, levels = c("PBS", "ASO 6"))
  
  par(mar=c(4,0,2,0))
  
  plot(
    NA, NA,
    xlim = c(0.5, 2.5),
    ylim = c(0, max(plot_elisa_rd$norm_conc) * 1.2),
    xlab = "", ylab = " ",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  axis(1, at = 1:2, labels = x_label[c("PBS", "ASO 6")], line=0,col = NA,  tick = FALSE )
  axis(1, at=0:3, labels=NA, tck=FALSE)
  
  axis(side = 2, at = 0:6/4, labels = FALSE, las = 2)
  
  stat <- plot_elisa_rd %>%
    group_by(treatment,col) %>%
    summarise(
      mean = mean(norm_conc, na.rm = TRUE),
      se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(treatment, levels = c("PBS", "ASO 6")))
    )
  
  
  
  with(stat, {
    suppressWarnings(arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col = col))
    segments(x - 0.3, mean, x + 0.3, mean, lwd = 3,col = col)
  })
  
  
  beeswarm(plot_elisa_rd$norm_conc ~ plot_elisa_rd$treatment,
           pch = 21, cex = 1.2,
           pwcol =plot_elisa_rd$col,
           add = TRUE)
  
  
  
  
  abline(h = 1, lty = 2, col = "gray50")
  abline(h = 0.5,  lty = 2, col = "gray80")
  
  
  
  legend(x='bottomleft', y=1, bty='n', cex=.8, lwd=2, legend=rev(unique(plot_elisa_rd$treatment)),
         col=rev(unique(plot_elisa_rd$col)),
         text.col=rev(unique(plot_elisa_rd$col)))
  
  mtext('CSF PrP', side = 3, line = 0.5, font = 2, cex = 0.8)
  
}





### Brain RNA ####




rna_rd = read_tsv('data/analytic/dyr_rna_ionis.tsv', col_types=cols())%>%
  filter(broad_study_id=='RD')%>%
  rename(expt=broad_study_id)%>%
  mutate(animal=as.character(animal))%>%
  select(-treatment)

plot_rna_rd <- rna_rd %>%
  left_join(blind %>% select(animal, cohort, genotype, inoculum,treatment), by=c('animal'))%>%
  mutate(col = case_when(
    genotype == "WT" & treatment=='ASO 6' ~ "#0001CD",
    genotype == "HET" ~ "#501378",
    genotype == "HOM" ~ "#A02422",
    treatment=='PBS' ~ "#A9A9A9"))



write_supp_table(plot_rna_rd%>%select(-col,-qpcr_percent,-cohort), 'RD Brain RNA individual')


aso_mean <- mean(plot_rna_rd$qpcr_normed[plot_rna_rd$treatment == "PBS"], na.rm = TRUE)



RD_Brain_RNA<- function(){ 
  par(xpd=F)
  
  n_by_group <- plot_rna_rd %>% count(treatment, name = "n")
  x_label <- setNames(
    paste0(levels(factor(plot_rna_rd$treatment, levels = c("PBS", "ASO 6"))),
           "\n(n=", n_by_group$n[match(c("PBS", "ASO 6"), n_by_group$treatment)], ")"),
    c("PBS", "ASO 6")
  )
  
  plot_rna_rd$treatment <- factor(plot_rna_rd$treatment, levels = c("PBS", "ASO 6"))
  
  
  par(mar=c(4,5,2,1))
  
  plot(
    NA, NA,
    xlim = c(0.5, 2.5),
    ylim = c(0, max(plot_rna_rd$qpcr_normed) * 1.3),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  mtext(side=2, line=3.2, text='residual (%)',cex=0.8)
  
  axis(1, at = 1:2, labels = x_label[c("PBS", "ASO 6")], line=0,col = NA,  tick = FALSE )
  axis(1, at=0:3, labels=NA, tck=FALSE)
  axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2)
  
  
  stat_rd <- plot_rna_rd %>%
    group_by(treatment,col) %>%
    summarise(
      mean = mean(qpcr_normed, na.rm = TRUE),
      se = sd(qpcr_normed, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(treatment, levels = c("PBS", "ASO 6")))
    )
  
  
  with(stat_rd, {
    arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col = col)
    segments(x - 0.3, mean, x + 0.3, mean, lwd = 3,col = col)
  })
  
  
  beeswarm(plot_rna_rd$qpcr_normed ~ plot_rna_rd$treatment,
           pch = 21, cex = 1.2,
           pwcol =plot_rna_rd$col,
           add = TRUE)
  
  abline(h = 1, lty = 2, col = "gray50")
  abline(h = 0.5,  lty = 2, col = "gray80")
  
  
  
  
  mtext('Brain RNA ', side = 3, line = 0.5, font = 2, cex = 0.8)
  
}

write_supp_table(plot_rna_rd %>%
                   group_by(treatment,col) %>%
                   summarise(
                     mean = mean(qpcr_normed, na.rm = TRUE),
                     se = sd(qpcr_normed, na.rm = TRUE) / sqrt(n()),
                     .groups = "drop"
                   ) %>% ungroup()%>%
                   mutate(
                     l95 = mean - 1.96 * se,
                     u95 = mean + 1.96 * se,
                     x = as.numeric(factor(treatment, levels = c("PBS", "ASO 6")))
                   )
                 %>%select(-col,-x), 'RD brain RNA individual')




## RG   ####

# all RG analyses exclude timed-sac animals and the N=1 incidental uninoculated animal (not part of original experimental design)
rg_keep_animals <- master %>%
  filter(expt == "RG", endpoint_comments != "timed sac") %>%
  pull(animal)

### survival ####
surv_ex <- survdata %>%
  filter(expt %in% c("RG"), animal %in% rg_keep_animals, cohort_disp != "Uninoculated")

write_supp_table(surv_ex %>% select(-color,-cohort), 'RG survival individual')


surv_smry<- surv_table(surv_ex, ref_group = "WT")%>%
  arrange(match(cohort_disp, c('WT','Het KO','Hom KO')))

write_supp_table(surv_smry, 'RG survival summarized')





RG_survival<- function(){ 
  par(mar=c(4,4,2,2))
  
  
  surv_ex <- survdata %>%
    filter(expt %in% c("RG"), animal %in% rg_keep_animals, cohort_disp != "Uninoculated") %>%
    mutate(
      cohort_disp = factor(
        cohort_disp,
        levels = c("WT", "Het KO", "Hom KO")
      )
    )
  
  
  sf = make_sf(surv_ex)
  
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
  
}




### weights ####




dat_ex <- weight_plot %>%
  filter(expt == "RG", animal %in% rg_keep_animals, cohort_disp != "Uninoculated") %>%
  select(cohort, animal, wtdpi, color, weight, peak, offpeak,cohort_disp) %>%
  ungroup()


write_supp_table(dat_ex%>%select(-color,-cohort), 'RG weights individual')

rel_weights_smry_RG<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(offpeak),
            sd = sd(offpeak),
            n = n(),
            l95 = mean(offpeak) - 1.96*sd(offpeak)/sqrt(n()),
            u95 = mean(offpeak) + 1.96*sd(offpeak)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry_RG%>%ungroup() %>%select(-color), 'RG weights summarized')



rg_weight_filtered <- weight_plot %>%
  filter(expt == 'RG', animal %in% rg_keep_animals, cohort_disp != "Uninoculated")

rg_peak_dpi_cohort <- rg_weight_filtered %>%
  group_by(cohort_disp, animal) %>%
  summarize(.groups = "drop", peak_dpi = wtdpi[which.max(weight)])

rg_p_vs_wt <- sapply(c("Het KO", "Hom KO"), function(g) {
  tryCatch(
    t.test(rg_peak_dpi_cohort$peak_dpi[rg_peak_dpi_cohort$cohort_disp == g],
           rg_peak_dpi_cohort$peak_dpi[rg_peak_dpi_cohort$cohort_disp == "WT"])$p.value,
    error = function(e) NA_real_
  )
})

weight_summary <- make_weight_summary(rg_weight_filtered, ref_group = "WT")%>%
  arrange(match(cohort_disp, c('WT','Het KO','Hom KO'))) %>%
  mutate(p_vs_wt = case_when(
    cohort_disp == "Het KO" ~ rg_p_vs_wt[["Het KO"]],
    cohort_disp == "Hom KO" ~ rg_p_vs_wt[["Hom KO"]],
    TRUE ~ NA_real_
  ))

write_supp_table(weight_summary, 'RG weights peak')


### peak weight dpi by animal and genotype ####

rg_peak_by_animal <- dat_ex %>%
  inner_join(blind %>% filter(expt == "RG") %>% distinct(animal, genotype), by = "animal") %>%
  arrange(animal, wtdpi) %>%
  group_by(cohort, cohort_disp, genotype, animal) %>%
  summarize(.groups = "drop",
            peak_weight = max(weight),
            peak_dpi    = wtdpi[weight == max(weight)][1])

write_supp_table(rg_peak_by_animal %>% select(-cohort), 'RG peak weight by animal')

rg_peak_dpi_genotype <- rg_peak_by_animal

p_vs_wt <- sapply(c("HE", "KO"), function(g) {
  tryCatch(
    t.test(rg_peak_dpi_genotype$peak_dpi[rg_peak_dpi_genotype$genotype == g],
           rg_peak_dpi_genotype$peak_dpi[rg_peak_dpi_genotype$genotype == "WT"])$p.value,
    error = function(e) NA_real_
  )
})

rg_peak_dpi_genotype_smry <- rg_peak_dpi_genotype %>%
  group_by(genotype) %>%
  summarize(.groups = "drop",
            n = n(),
            mean_peak_dpi = mean(peak_dpi),
            l95 = mean(peak_dpi) - 1.96 * sd(peak_dpi) / sqrt(n()),
            u95 = mean(peak_dpi) + 1.96 * sd(peak_dpi) / sqrt(n())) %>%
  mutate(p_vs_wt = case_when(
    genotype == "HE" ~ p_vs_wt[["HE"]],
    genotype == "KO" ~ p_vs_wt[["KO"]],
    TRUE ~ NA_real_
  )) %>%
  arrange(match(genotype, c("WT", "HE", "KO")))

write_supp_table(rg_peak_dpi_genotype_smry, 'RG peak weight dpi by genotype')


RG_weight<- function(){
  
  par(mar=c(3,4,2,2))
  plot_weights_max(rel_weights_smry_RG,ylims=c(-1,0.2))
  
  
  
  
}







## RE  ####
### survival ####


surv_ex <- survdata %>%
  filter(expt %in% c("RE"))%>%
  filter(!grepl("-died", cohort))%>%
  mutate(cohort_disp=factor(cohort_disp, levels=c("Saline","ASO 6 1000ug")))


write_supp_table(surv_ex %>% select(-color,-cohort), 'RE survival individual')

surv_smry<- surv_table(surv_ex, ref_group = "Saline")

write_supp_table(surv_smry, 'RE survival summarized')



RE_survival<-function(){
  
  par(mar=c(3,4,2,1))
  
  
  surv_ex <- survdata %>%
    filter(expt %in% c("RE"))%>%
    filter(!grepl("-died", cohort))%>%
    mutate(cohort_disp=factor(cohort_disp, levels=c("Saline","ASO 6 1000ug")))
  
  anno_date <- as.numeric(unique(na.omit(surv_ex$txdpi)))
  
  sf = make_sf(surv_ex)
  plot_sf(sf,xlims = c(0, 350))
  
  points(
    x = anno_date,
    y = 1.04,
    pch = 25,                     
    bg = 'black',        
    col = 'black',      
    cex = 0.5
  )
  
  
}

## RP  ####


rat_it_icv %>%
  group_by(region, route, dose_ug) %>%
  summarize(.groups='keep',
            n = n(),
            mean = mean(residual),
            l95 = lower(residual),
            u95 = upper(residual)) %>%
  ungroup() -> rat_it_icv_smry

write_supp_table(rat_it_icv_smry, 'Summary of ASO 6 potency in rat by IT catheter vs. ICV route of administration.')


### survival ####
surv_ex <- survdata %>%
  filter(expt %in% c("RP"))%>%
  filter(!grepl("-died", cohort))%>%
  mutate(cohort_disp=factor(cohort_disp, levels=c("Vehicle","ASO 6 300ug", "ASO 6 600ug"  ,"ASO 6 1000ug")))


write_supp_table(surv_ex %>% select(-color,-cohort), 'RP survival individual')

surv_smry<- surv_table(surv_ex, ref_group = "Vehicle")

write_supp_table(surv_smry, 'RP survival summarized')

# log-rank test, RP ASO 6 600 ug vs. 1000 ug (not reported in Tables S1/S2)
rp_600_vs_1000 <- survdiff(Surv(dpi, acm) ~ cohort_disp,
                           data = surv_ex %>% filter(cohort_disp %in% c('ASO 6 600ug', 'ASO 6 1000ug')) %>%
                             mutate(cohort_disp = droplevels(factor(cohort_disp))))
rp_600_vs_1000_p <- 1 - pchisq(rp_600_vs_1000$chisq, df = 1)

RP_survival<- function(){ 
  
  
  par(mar=c(4,4,2,2))
  
  
  surv_ex <- survdata %>%
    filter(expt %in% c("RP"))%>%
    filter(!grepl("-died", cohort))%>%
    mutate(cohort_disp=factor(cohort_disp, levels=c("Vehicle","ASO 6 300ug", "ASO 6 600ug"  ,"ASO 6 1000ug")))
  
  sf = make_sf(surv_ex)
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi) + 10))
  
  
}




### weights ####



dat_ex<- weight_plot %>%
  filter(expt == 'RP') %>%
  filter(!grepl("-died", cohort))%>%
  select(cohort, animal, wtdpi, color, weight, peak, offpeak,cohort_disp) %>%
  ungroup()%>%
  mutate(cohort_disp=factor(cohort_disp, levels=c("Vehicle","ASO 6 300ug", "ASO 6 600ug"  ,"ASO 6 1000ug")))

write_supp_table(dat_ex%>%select(-color,-cohort), 'RP weights peak')

rel_weights_smry_RP<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(offpeak),
            sd = sd(offpeak),
            n = n(),
            l95 = mean(offpeak) - 1.96*sd(offpeak)/sqrt(n()),
            u95 = mean(offpeak) + 1.96*sd(offpeak)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry_RP%>%ungroup() %>%select(-color), 'RP weights summarized')



weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt == 'RP'), ref_group = "Vehicle")%>%
  arrange(match(cohort_disp, c('Vehicle','ASO 6 300ug','ASO 6 600ug','ASO 6 1000ug')))

write_supp_table(weight_summary, 'RP cohort weight summary')




RP_weight<-function(){
  par(mar=c(3,4,2,2))
  
  
  plot_weights_max(rel_weights_smry_RP,xlims = c(0,330),ylims=c(-1,0.2))
  
}






# HAMSTER ####

## ELISA ####


hamster_elisa <- function(){ 
  par(mar = c(3,4,2,2))
  
  geno_labels = c("WT","Het KO","Hom KO")
  leg<- plot_csELISA(ana_132, "Hamster", cex_a=0.65, geno_labels=geno_labels)
  
  legend(x='topright', y=1, bty='n', cex=.5, lwd=2, legend=geno_labels[match(leg$genotype, c("WT","HET","HOM"))],
         col=leg$col,
         text.col=leg$col)
  
  
  
  mtext(side=2, line=3.1, text='residual PrP (%)',cex=0.8)
  mtext('brain', side = 3, line = 0.5, font = 1, cex = 0.8)
  
}

stat_Hamster <- {
  df <- ana_132 %>%
    filter(genotype %in% c("WT", "HET", "HOM")) %>%
    mutate(col = case_when(genotype == "WT"  ~ "#0001CD",
                           genotype == "HET" ~ "#501378",
                           genotype == "HOM" ~ "#A02422"))
  wt_mean <- mean(df$ngml_av[df$genotype == "WT"], na.rm = TRUE)
  df %>%
    mutate(norm_conc = ngml_av / wt_mean) %>%
    group_by(genotype, col) %>%
    summarise(mean = mean(norm_conc, na.rm = TRUE),
              se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    ungroup() %>%
    mutate(l95 = mean - 1.96 * se,
           u95 = mean + 1.96 * se,
           x = as.numeric(factor(genotype, levels = c("WT","HET","HOM"))))
}
write_supp_table(stat_Hamster %>% select(-x,-col), "Hamster ELISA summarized")


##  ASO ####


elisa_hamster = elisas %>%
  filter(is.na(exclude),
         species=='hamster')%>%
  mutate(animal_id=as.character(unlist(animal_id)),
         sample=unlist(sample),
         ngml_av=unlist(ngml_av),
         se_mean=unlist(se_mean))%>%
  rename(animal=animal_id,expt=experiment)

plot_elisa_hamster <- elisa_hamster %>%
  left_join(blind %>% select(animal, cohort, genotype, treatment, inoculum), by=c('animal'))%>%
  mutate(col = case_when(
    genotype == "WT" & treatment=='ASO 1' ~ "#CFB53B",
    genotype == "HET" ~ "#501378",
    genotype == "HOM" ~ "#A02422",
    treatment=='none' ~ "#A9A9A9"))

aso_mean <- mean(plot_elisa_rd$ngml_av[plot_elisa_rd$treatment == "PBS"], na.rm = TRUE)

plot_elisa_hamster <- plot_elisa_hamster %>%
  group_by(tissue) %>%
  mutate(
    none_mean = mean(ngml_av[treatment == "none"], na.rm = TRUE),
    norm_conc = ngml_av / none_mean 
  ) %>%
  ungroup()

write_supp_table(plot_elisa_hamster%>% select(-col,-exclude,-cohort,-none_mean), 'Hamster ASO individual')


plot_elisa_hamster_smry<- plot_elisa_hamster %>%
  group_by(tissue, treatment) %>%  
  summarise(
    N = n(),
    mean_ngml = mean(ngml_av, na.rm = TRUE),
    sd_ngml   = sd(ngml_av, na.rm = TRUE),
    se_ngml   = sd(ngml_av, na.rm = TRUE) / sqrt(N),
    mean_norm = mean(norm_conc, na.rm = TRUE),
    sd_norm   = sd(norm_conc, na.rm = TRUE),
    se_norm   = sd(norm_conc, na.rm = TRUE) / sqrt(N),
    .groups = "drop"
  ) %>% ungroup()%>%
  mutate(
    mean_ngml = round(mean_ngml, 2),
    sd_ngml   = round(sd_ngml, 2),
    se_ngml   = round(se_ngml, 3),
    mean_norm = round(mean_norm, 2),
    sd_norm   = round(sd_norm, 2),
    se_norm   = round(se_norm, 3)
  )%>%
  arrange(match(treatment, c('none','ASO 1')))%>%
  arrange(match(tissue, c('BH','CSF')))

write_supp_table(plot_elisa_hamster_smry, 'Hamster ASO summarized')



hamster_aso<- function(tis){
  par(mar = c(3,4,2,2))
  
  if(tis=='BH') {
    title = 'brain'
    used_dilution = 200 # Meredith ran 1:20 of a 10% wt/vol
  } else if (tis=='CSF') {
    title = 'CSF'
    used_dilution = 15 # Meredith ran 1:15 and 1:40, so 1:15 defines the LLQ
  }
  
  
  df<- plot_elisa_hamster %>%
    filter(tissue==tis)
  
  n_by_group <- df %>%count(treatment, name = "n")
  x_label <- setNames(
    paste0(levels(factor(df$treatment, levels = c("none", "ASO 1"))),
           "\n(n=", n_by_group$n[match(c("none", "ASO 1"), n_by_group$treatment)], ")"),
    c("none", "ASO 1")
  )
  
  df$treatment <- factor(df$treatment, levels = c("none", "ASO 1"))
  
  
  
  
  
  plot(
    NA, NA,
    xlim = c(0.5, 2.5),
    ylim = c(0, max(df$norm_conc) * 1.2),
    xlab = "", ylab = "residual PrP (%)",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  axis(1, at = 1:2, labels = x_label[c("none", "ASO 1")], line=0,col = NA,  tick = FALSE ,cex.axis=0.8)
  axis(1, at=0:3, labels=NA, tck=FALSE)
  
  axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2)
  
  untreated_mean = mean(df$ngml_av[df$treatment=='none'])
  LLQ = round(0.05 * used_dilution / untreated_mean, 5)
  stat <- df %>%
    group_by(treatment,col) %>%
    summarise(
      mean = mean(norm_conc, na.rm = TRUE),
      se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(treatment, levels = c("none", "ASO 1")))
    )
  
  
  with(stat, {
    arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col = col)
    segments(x - 0.25, mean, x + 0.25, mean, lwd = 3,col = col)
  })
  
  
  beeswarm(df$norm_conc ~ df$treatment,
           pch = 21, cex = 1.2,
           pwcol =df$col,
           add = TRUE)
  
  
  
  
  abline(h = 1, lty = 2, col = "gray50")
  abline(h = .5,  lty = 2, col = "gray80")
  abline(h = LLQ, lty = 2, col = "red")
  mtext(side=4, at=LLQ, line=0.15, las=2, text='LLQ', col='red', cex=0.7)
  
  yleg = ifelse(tis=='CSF',1.6,1.5)
  par(xpd=T)
  legend(x=1.7, y=yleg, bty='n', cex=.65, lwd=2, legend=unique(df$treatment),
         col=unique(df$col),
         text.col=unique(df$col))
  par(xpd=F)
  mtext(title, side = 3, line = 0.5, font = 1, cex = 0.8)
  

  
  
  
}





## SG/SH ####

### survival ####
surv_ex <- survdata %>%
  filter(expt %in% c("SG", "SH")) 


write_supp_table(surv_ex %>% select(-color,-cohort), 'SG/SH survival individuals')



surv_smry<- surv_table(surv_ex, ref_group = "Untreated WT 263K")%>%
  arrange(match(cohort_disp, c('Untreated WT 263K','WT uninoculated','Hom KO 263K','Het KO 263K')))

write_supp_table(surv_smry, 'SG/SH survival summarized')


SG_SH_survival<-function(){
  par(mar=c(4,4,2,2))
  
  surv_ex <- survdata %>%
    filter(expt %in% c("SG", "SH")) 
  
  
  sf = make_sf(surv_ex) 
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
  
}


###weight ####

dat_ex<- weight_plot %>%
  filter(expt %in% c("SG", "SH")) %>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  


write_supp_table(dat_ex%>%select(-color,-cohort) , 'SG/SH weights individual')


rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%
  select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'SG/SH weights summarized')


weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt %in% c("SG", "SH"))  , ref_group = "Untreated WT 263K")%>%
  arrange(match(cohort_disp, c('Untreated WT 263K','WT uninoculated','Hom KO 263K','Het KO 263K')))

write_supp_table(weight_summary, 'SG/SH weights peak')







## SC/SD/SE/SK  ####


### survival ####
surv_ex <- survdata %>%
  filter(expt %in% c('SC','SD','SE','SK'),
         cohort %in% c('Hamster-HE-Untreated', 'Hamster-WT-Untreated', 'Hamster-WT-ASO1-300ug-64dpi','Hamster-WT-ASO1-300ug-43dpi')) %>%
  filter(!is.na(cohort))



write_supp_table(surv_ex %>% select(-color,-cohort), 'SC/SD/SE/SK survival individual')


sf = make_sf(surv_ex)

surv_smry<- surv_table(surv_ex, ref_group = "Untreated WT 263K")%>%
  arrange(match(cohort_disp, c('Untreated WT 263K','Het KO 263K','WT 263K ASO 1 at 43 dpi','WT 263K ASO 1 at 64 dpi')))

write_supp_table(surv_smry, 'SC/SD/SE/SK survival summarized')


# Pooled log-rank test reproducing P=0.029 (Figure 2C): ASO 1 treatment
# (either dosing timepoint) vs untreated, restricted to WT 263K hamsters.
sdobj = survdiff(Surv(dpi, acm) ~ treatment, data = surv_ex %>% filter(genotype == 'WT'))
sdobj_tbl <- tibble(
  group    = sub('.*=', '', names(sdobj$n)),
  N        = as.integer(sdobj$n),
  observed = as.numeric(sdobj$obs),
  expected = as.numeric(sdobj$exp),
  chisq    = as.numeric(sdobj$chisq),
  df       = length(sdobj$n) - 1L,
  p_value  = 1 - pchisq(sdobj$chisq, length(sdobj$n) - 1L)
)
write_supp_table(sdobj_tbl, 'Log-rank test of ASO 1 treatment vs untreated in WT 263K hamsters')




SC_SD_SE_SK_survival<- function(){
  par(mar=c(4,4,2,2))
  surv_ex <- survdata %>%
    filter(expt %in% c('SC','SD','SE','SK'),
           cohort %in% c('Hamster-HE-Untreated', 'Hamster-WT-Untreated', 'Hamster-WT-ASO1-300ug-64dpi','Hamster-WT-ASO1-300ug-43dpi')) %>%
    filter(!is.na(cohort))
  sf = make_sf(surv_ex)
  
  anno_dpi= surv_ex %>% filter(!is.na(txdpi)) %>% distinct(txdpi) %>% pull() %>% as.integer()
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)+10),anno_dpi= anno_dpi)
  
  sdobj = survdiff(Surv(dpi, acm) ~ treatment, data=surv_ex %>%
                     filter(genotype=='WT'))
  
}



### weight #### 


dat_ex<- weight_plot %>%
  filter(expt %in% c('SC','SD','SE','SK'),
         cohort %in% c('Hamster-HE-Untreated', 'Hamster-WT-Untreated', 'Hamster-WT-ASO1-300ug-64dpi','Hamster-WT-ASO1-300ug-43dpi')) %>%
  filter(!is.na(cohort))%>%
  select(cohort, animal, wtdpi, color, weight, peak, cohort_disp,offpeak) %>%
  ungroup()

# Figure S3C normalizes to each animal's weight at ~50 dpi (the weighing nearest
# 50 dpi, within 5 days; animals with none are dropped) rather than to peak weight.
dat_ex <- dat_ex %>%
  group_by(animal) %>%
  mutate(baseline50_dpi = wtdpi[which.min(abs(wtdpi - 50))],
         baseline50_dpi = ifelse(abs(baseline50_dpi - 50) <= 5, baseline50_dpi, NA_integer_),
         baseline50 = weight[match(baseline50_dpi[1], wtdpi)],
         rel50 = (weight - baseline50) / baseline50) %>%
  ungroup()

write_supp_table(dat_ex%>%select(-color,-cohort,-baseline50_dpi) , 'SC/SD/SE/SK weights individual')





rel_weights_smry_SC_SD_SE_SK <- dat_ex %>%
  filter(!is.na(rel50), wtdpi >= baseline50_dpi) %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel50),
            sd = sd(rel50),
            n = n(),
            l95 = mean(rel50) - 1.96*sd(rel50)/sqrt(n()),
            u95 = mean(rel50) + 1.96*sd(rel50)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)


write_supp_table(rel_weights_smry_SC_SD_SE_SK%>%ungroup() %>%select(-color), 'SC/SD/SE/SK weights summarized')




weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt %in% c('SC','SD','SE','SK'),     
                                               cohort %in% c('Hamster-HE-Untreated', 
                                                             'Hamster-WT-Untreated', 
                                                             'Hamster-WT-ASO1-300ug-64dpi',
                                                             'Hamster-WT-ASO1-300ug-43dpi')), 
                                      ref_group = "Untreated WT 263K")

write_supp_table(weight_summary%>%
                   arrange(match(cohort_disp, c('Untreated WT','Untreated HE','WT 43dpi','WT 63dpi'))), 'SC/SD/SE/SK weights peak')






SC_SD_SE_SK_weight<-function(){
  
  par(mar=c(4,4,2,2))
  
  plot_weights_max(rel_weights_smry_SC_SD_SE_SK,xlims=c(50, 200),ylims = c(-0.4,0.8), euth_level=-0.15, ylab='weight vs. 50 dpi', legend_pos='topleft')
  axis(side=1, at=50, lwd=0, labels=50, line=-0.5)
  
}





### nest ####
dat_ex <- nest_plot %>%
  filter(expt %in% c('SC','SD','SE','SK'),
         cohort %in% c('Hamster-HE-Untreated', 'Hamster-WT-Untreated', 'Hamster-WT-ASO1-300ug-64dpi','Hamster-WT-ASO1-300ug-43dpi')) %>%
  filter(!is.na(cohort))

write_supp_table(dat_ex%>%select(-color,-cohort), 'SC/SD/SE/SK nests individual')

SC_SD_SE_SK_nests<-function(){
  
  par(mar=c(4,5,2,4))
  
  anno <- blind %>% filter(expt %in% c('SC','SD','SE','SK')) %>% pull(treatment_start_dpi) %>% unique()
  anno <- anno[!is.na(anno)]
  
  
  
  plot_nest(dat_ex,
            anno_dpi = if (length(anno)) anno else NULL,
            xlims = c(0, 200))
  
  mtext('SC/SD/SE/SK', side = 3, line = 0.5, font = 2, cex = 0.9)
}



## SH/SB/SC/SD/SE/SK RT-QuIC  ####



hamster_rtq = read_tsv('data/analytic/dyr_rtquics.tsv', col_types=cols()) %>%
  filter(experiment %in% c('SH','SB','SC','SD','SE','SK'),
         is.na(exclude)) %>%
  clean_names() %>% 
  mutate(animal_id=as.character(unlist(animal_id)) ,
         sample=unlist(sample))%>%
  rename(animal=animal_id,expt=experiment) %>%
  arrange(tissue, desc(rta), animal, dilution) %>%
  group_by(animal, tissue, dilution) %>% slice(1) %>% ungroup() %>% # for samples re-run, remove the earlier version
  arrange(rta, animal, tissue, desc(dilution)) %>%
  inner_join(master %>% select(animal,dob,inoculation_date,death_date), by=c('animal'))%>%
  left_join(blind %>% select(cohort,animal,expt), by=c('animal','expt'))%>%
  left_join(meta,by=c('cohort','expt','species'))%>%
  mutate(days=as.integer(death_date-inoculation_date))%>%
  filter( !is.na(cohort_disp))



rtq_plot_hamster <- hamster_rtq %>%
  filter(tissue=='brain') %>%
  group_by(sample, days, dilution, color, cohort_disp, expt) %>%
  summarize(
    positive_replicates = sum(positive_replicates),
    total_replicates    = sum(total_replicates),
    .groups = "drop"
  )%>% ungroup() %>%
  arrange(sample, days, desc(dilution)) %>%
  group_by(sample, days,color,cohort_disp,expt) %>%
  summarize(
    sk = spearman_karber(
      x = positive_replicates,
      n = total_replicates,
      d = log10(dilution)
    ),
    .groups = "drop"
  )%>% ungroup() %>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2], llq=sk[,3] + -log10(bh_unit_factor), ulq=sk[,4] + -log10(bh_unit_factor)) %>% 
  select(-sk) %>%
  ungroup() 



brain_hamster_titer_model = lm(titer ~ days : cohort_disp, data=rtq_plot_hamster)
summary(brain_hamster_titer_model)$coefficients %>%
  as_tibble(rownames='variable') %>%
  clean_names() %>%
  rename(p = pr_t) %>%
  mutate(daily_change = 10^estimate-1) %>%
  mutate(doubling_time = log10(2) / estimate) -> hamster_brain_titer_model_coefficients

brain_hamster_titer_model_wt_het_only = lm(titer ~ days : cohort_disp, data=rtq_plot_hamster %>% filter(cohort_disp %in% c('Untreated WT 263K','Het KO 263K')))

write_supp_table(hamster_brain_titer_model_coefficients, 'Linear model of hamster brain RT-QuIC titer.')

write_supp_table(rtq_plot_hamster%>%select(-color), 'Hamster RT-quic individual')


rtq_expt_counts <- rtq_plot_hamster %>%
  mutate(approx_days = round(days / 15) * 15) %>%
  group_by(approx_days, cohort_disp, expt) %>%
  summarise(
    n_expt = n(), 
    .groups = "drop"
  ) %>% ungroup()%>%
  tidyr::pivot_wider(
    names_from = expt,
    values_from = n_expt,
    names_prefix = "n_expt_"
  )

hamster_rtq_smry <- rtq_plot_hamster %>%
  select(-color) %>%
  mutate(approx_days = round(days / 15) * 15) %>%
  group_by(approx_days, cohort_disp) %>%
  summarise(
    n     = n(),
    mean_titer = mean(titer, na.rm = TRUE),
    se_titer   = sd(titer, na.rm = TRUE) / sqrt(n),
    l95   = mean_titer - 1.96 * se_titer,
    u95   = mean_titer + 1.96 * se_titer,
    .groups = "drop"
  ) %>% ungroup()%>%
  left_join(rtq_expt_counts, by = c("approx_days", "cohort_disp"))


write_supp_table(hamster_rtq_smry, 'Hamster RT-quic summarized by days')


rtq_plot_hamster_csf <- hamster_rtq %>%
  filter(tissue=='CSF') %>%
  group_by(sample, days, dilution,color,cohort_disp,expt) %>%
  summarize(
    positive_replicates = sum(positive_replicates),
    total_replicates    = sum(total_replicates),
    .groups = "drop"
  )%>% ungroup() %>%
  arrange(sample, days, desc(dilution)) %>%
  group_by(sample, days,color,cohort_disp,expt) %>%
  summarize(
    sk = spearman_karber(
      x = positive_replicates,
      n = total_replicates,
      d = log10(dilution)
    ),
    .groups = "drop"
  )%>% ungroup() %>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2], llq=sk[,3] + -log10(bh_unit_factor), ulq=sk[,4] + -log10(bh_unit_factor)) %>% 
  select(-sk) %>%
  ungroup() 



csf_hamster_titer_model = lm(titer ~ days : cohort_disp, data=rtq_plot_hamster_csf)
summary(csf_hamster_titer_model)$coefficients %>%
  as_tibble(rownames='variable') %>%
  clean_names() %>%
  rename(p = pr_t) -> hamster_csf_titer_model_coefficients

write_supp_table(hamster_csf_titer_model_coefficients, 'Linear model of hamster CSF RT-QuIC titer.')


rtq_expt_counts_csf <- rtq_plot_hamster_csf %>%
  mutate(approx_days = round(days / 15) * 15) %>%
  group_by(approx_days, cohort_disp, expt) %>%
  summarise(
    n_expt = n(), 
    .groups = "drop"
  ) %>% ungroup()%>%
  tidyr::pivot_wider(
    names_from = expt,
    values_from = n_expt,
    names_prefix = "n_expt_"
  )

hamster_rtq_smry_csf <- rtq_plot_hamster_csf %>%
  select(-color) %>%
  mutate(approx_days = round(days / 15) * 15) %>%
  group_by(approx_days, cohort_disp) %>%
  summarise(
    n     = n(),
    mean_titer = mean(titer, na.rm = TRUE),
    se_titer   = sd(titer, na.rm = TRUE) / sqrt(n),
    l95   = mean_titer - 1.96 * se_titer,
    u95   = mean_titer + 1.96 * se_titer,
    .groups = "drop"
  ) %>% ungroup()%>%
  left_join(rtq_expt_counts, by = c("approx_days", "cohort_disp"))


SC_SD_SE_SK_rtquic<- function(tissue='brain') {
  
  par(mar=c(4,4,2,2))
  
  if (tissue=='brain') {
    plot_data = rtq_plot_hamster
    ylims = c(5,13)
    ylab = 'titer (log10 seeds/g)'
  } else if (tissue=='CSF') {
    plot_data = rtq_plot_hamster_csf
    ylims = c(2.5,8)
    ylab = 'titer (log10 seeds/mL)'
  }
  
  
  
  
  par(xpd=F)
  xlims = c(0,180)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
  
  axis(side=1, at=seq(0,200,30), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,200,30), lwd=0, labels=seq(0,200,30), line=-0.5)
  axis(side=1, at=seq(0,200,30), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='dpi',cex=0.8)
  yats=log10(rep(1:9,12) * rep(10^(1:12),each=9))
  ybigs=1:13
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=3, cex=0.8, text=tissue)
  
  mtext(side=2, line=2.75, text=ylab,cex=0.8)
  
  par(xpd=T)
  points(plot_data$days,plot_data$titer,col=plot_data$color)
  par(xpd=F)
  
  for (g in unique(plot_data$cohort_disp)) {
    
    subs<- plot_data%>%filter(cohort_disp==g)
    lm_model = lm(titer ~ days, data=subs)
    
    x = seq(0:180)
    y <- predict(lm_model, data.frame(days = x))
    
    points(x, y, type='l', col=subs$color[1], lwd=2)
  }
  par(xpd=F)
  

  abline(h=min(plot_data$llq), lty=3, col='red')
  mtext(side=4, at=min(plot_data$llq), col='red', las=2, text='LLQ',cex=0.7)
  
  
}



## SC_SD_SE_SK NfL  ####




nfl_SC_SD_SE_SK = ellas %>%
  filter(experiment %in% c('SH','SB','SC','SD','SE','SK','SG'),
         is.na(exclude)) %>%
  filter(gn_rs == 3) %>%
  rename(animal=animal_id,expt=experiment)%>%
  left_join(master %>% select(expt,animal,dob,death_date,inoculation_date), by=c('expt','animal'))%>%
  mutate(sample_date=as.Date(sample_date),
         sample_dpi = if_else(
           is.na(sample_dpi),
           as.numeric(sample_date - inoculation_date),
           sample_dpi)) %>%
  filter(expt != "SH" | sample_date != death_date)%>%
  left_join(blind %>% select(expt,cohort,animal,genotype), by=c('expt','animal'))%>%
  left_join(meta %>% distinct(expt,cohort,color,cohort_disp), by=c('expt','cohort'))%>%
  select(animal,species,genotype,expt,inoculation_date,tissue,death_date,sample_date,sample_dpi,sample_name,mean_rfu,rfu_percent_cv,mean_conc,conc_percent_cv,color,cohort_disp)%>%
  mutate(
    approx_sample_dpi = round(sample_dpi / 15) * 15
  ) %>%
  mutate(mean_conc = as.numeric(mean_conc))






plot_nfl_SC_SD_SE_SK<- nfl_SC_SD_SE_SK %>%
  filter(!is.na(animal))%>%
  group_by(cohort_disp,approx_sample_dpi,color) %>%
  summarise(
    mean = mean(mean_conc, na.rm = TRUE),
    sd   = sd(mean_conc, na.rm = TRUE),
    n = n(),
    l95 = mean- 1.96*sd(mean_conc)/sqrt(n()),
    u95 = mean + 1.96*sd(mean_conc)/sqrt(n()),
    .groups   = "drop"
  ) %>% ungroup()%>%
  filter(!is.na(cohort_disp))%>%
  select(cohort_disp,approx_sample_dpi,everything())


write_supp_table(plot_nfl_SC_SD_SE_SK%>%select(-color), 'SC/SD/SE/SK NfL summarized')




ks_by_dpi <- nfl_SC_SD_SE_SK %>%
  filter(cohort_disp %in% c("Untreated WT 263K", "WT uninoculated")) %>%
  group_by(approx_sample_dpi) %>%
  summarize(
    n_263K   = sum(cohort_disp == "Untreated WT 263K"),
    n_uninoc = sum(cohort_disp == "WT uninoculated"),
    p_value = if (n_263K >= 2 & n_uninoc >= 2) {
      tryCatch(
        ks.test(
          mean_conc[cohort_disp == "Untreated WT 263K"],
          mean_conc[cohort_disp == "WT uninoculated"]
        )$p.value,
        error = function(e) NA
      )
    } else {
      NA
    },
    .groups = "drop"
  )%>% ungroup()%>%
  mutate(
    sig = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )



write_supp_table(ks_by_dpi%>%select(-sig), 'Hamster NfL ks test')


SC_SD_SE_SK_nfl <- function(){ 
  par(mar=c(4,4,2,2))
  
  xlims <-c(0,170)
  ylims =c(0,2000)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
  
  
  axis(side=1, at=seq(0,200,30), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,200,30), lwd=0, labels=seq(0,200,30), line=-0.5)
  axis(side=1, at=seq(0,200,30), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=seq(0,2000,500), lwd=0, labels=seq(0,2000,500), line=-0.5, las=2)
  axis(side=2, at=seq(0,2000,100), labels=NA, tck=-0.025)
  
  mtext(side=2, line=2.75, text='NfL (pg/ml)',cex=0.8)
  par(xpd=T)
  points(nfl_SC_SD_SE_SK$sample_dpi,nfl_SC_SD_SE_SK$mean_conc,col=nfl_SC_SD_SE_SK$color)
  par(xpd=F)
  
  for (co in unique(plot_nfl_SC_SD_SE_SK$cohort_disp)) {
    subs <- plot_nfl_SC_SD_SE_SK %>%
      arrange(approx_sample_dpi)%>%
      filter(cohort_disp==co,
             n>1)
    
    
    
    polygon(x=c(subs$approx_sample_dpi,rev(subs$approx_sample_dpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,0.5),border=NA)
    lines(subs$approx_sample_dpi, subs$mean, lwd = 2.5, col = subs$color)
    points(subs$approx_sample_dpi, subs$mean, pch = 16, col = subs$color)
    
  }
  
  
  star_pos <- plot_nfl_SC_SD_SE_SK %>%
    group_by(approx_sample_dpi) %>%
    summarize(
      y = max(u95, mean, na.rm = TRUE),
      .groups = "drop"
    )%>% ungroup()
  
  
  anno <- left_join(ks_by_dpi, star_pos, by = "approx_sample_dpi")
  
  
  text(
    x = anno$approx_sample_dpi,
    y = anno$y +250,  
    labels = anno$sig,
    cex = 1,
    font = 2
  )
  
  plot_nfl_SC_SD_SE_SK  %>%
    distinct(cohort_disp, color) %>%
    arrange(color) -> leg
  
  legend(1,2150, bty='n', cex=.7, lwd=2, legend=leg$cohort_disp,
         col=leg$color,
         text.col=leg$color)
  
  mtext('plasma NfL', side = 3, line = 0.5, font = 1, cex = 0.8)
  
}










## SZ  ####

### RT_quic ####
rt_quic_sz = read_tsv('data/analytic/dyr_rtquics.tsv', col_types=cols()) %>%
  filter(experiment=='SZ',
         is.na(exclude))%>%
  clean_names()%>% 
  mutate(animal_id=as.character(unlist(animal_id)),
         sample=unlist(sample))%>%
  rename(animal=animal_id,expt=experiment)


plot_rt_quic <-  rt_quic_sz %>%
  left_join(blind %>% select(animal, cohort, genotype, treatment, inoculum), by=c('animal'))%>%
  left_join(master %>%
              select(-cage) %>%
              mutate(dpi = as.integer(death_date - inoculation_date)), by=c('animal','expt'))%>%
  left_join(gt_meta , by='genotype')%>%
  select(expt, species, animal, tissue,inoculum, genotype, color, dpi, treatment, dilution, positive_replicates, total_replicates)%>%
  mutate(dpi=as.factor(dpi),
         group=case_when(
           dpi == "43" &treatment=='none'~  '43dpi none',
           dpi == "73" &treatment=='none'~  '73dpi none',
           dpi == "73" &treatment=='ASO 1'~  '73dpi ASO 1',
         ),
         color = case_when(
           group == "43dpi none" ~ "#0001CD",
           group == "73dpi none" ~ "#fb8500",
           group == "73dpi ASO 1" ~ "#CFB53B",
           TRUE ~ color))%>%
  filter(!is.na(group))


for (tis in unique(plot_rt_quic$tissue)) {
  
  df <- plot_rt_quic %>% filter(tissue == tis)
  df$group <- factor(df$group, levels = c("43dpi none", "73dpi none", "73dpi ASO 1"))
  
  df %>%
    group_by(expt, group, species, animal, inoculum, genotype, dpi, treatment, color) %>%
    summarize(.groups='keep',
              sk = spearman_karber(x = positive_replicates,
                                   n = total_replicates,
                                   d = log10(dilution))) %>% ungroup()%>%
    mutate(titer = sk[,1] + -log10(bh_unit_factor), se = sk[,2]) %>%
    select(-sk) %>% ungroup() -> titers
  
  write_supp_table(titers %>% select(-color),
                   paste0("SZ ", tis, " RT-quic individual"))
  
  stat <- titers %>%
    group_by(group, color) %>%
    summarise(
      n = n(),
      mean = mean(titer, na.rm = TRUE),
      se = sd(titer, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      l95 = mean - 1.96 * se,
      u95 = mean + 1.96 * se,
      x = as.numeric(factor(group))
    )
  
  write_supp_table(stat %>% select(-color, -x),
                   paste0("SZ ", tis, " RT-quic summarized"))
  t_tests <- tibble::tibble(
    comparison = c("73dpi none vs. 43dpinone", "73dpi ASO 1 vs 43dpi none", "73dpi ASO 1 vs 73dpi none"),
    p_value = c(
      tryCatch(
        t.test(titers$titer[titers$group=="73dpi none"],
               titers$titer[titers$group=="43dpi none"])$p.value,
        error = function(e) NA_real_
      ),
      tryCatch(
        t.test(titers$titer[titers$group=="73dpi ASO 1"],
               titers$titer[titers$group=="43dpi none"])$p.value,
        error = function(e) NA_real_
      ),
      tryCatch(
        t.test(titers$titer[titers$group=="73dpi ASO 1"],
               titers$titer[titers$group=="73dpi none"])$p.value,
        error = function(e) NA_real_
      )
    )
  )
  
  write_supp_table(t_tests,
                   paste0("SZ ", tis, " RT-quic T-tests"))
}



SZ_rtqic<- function(){
  
  for (tis in unique(plot_rt_quic$tissue)) {
    
    df<- plot_rt_quic %>%
      filter(tissue==tis)
    
    
    df$group <- factor(df$group, levels = c("43dpi none", "73dpi none","73dpi ASO 1"))
    
    
    
    
    if (tis == "brain") {
      par(mar=c(4,4,2,1))
      
      ylims=c(9,12)
      df %>%
        group_by(expt, group,species, animal, inoculum, genotype, dpi, treatment, color) %>%
        summarize(.groups='keep',
                  sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution))) %>% ungroup()%>%
        mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2]) %>%
        select(-sk) %>% ungroup() -> titers
      
      
    }else { ### tissue == CSF
      
      par(mar=c(4,1,2,2))
      
      ylims=c(4,11)
      
      df %>%
        group_by(expt, group,species, animal, inoculum, genotype, dpi, treatment, color) %>%
        summarize(.groups='keep',
                  sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution))) %>% ungroup()%>%
        mutate(titer = sk[,1] + -log10(0.001), se=sk[,2]) %>%
        select(-sk) %>% ungroup() -> titers
      
    }
    
    
    
    
    
    plot(
      NA, NA,
      xlim = c(0.5, 3.5),
      ylim = ylims,
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n",
      xaxs = "i", yaxs = "i",
      bty  = "n"
    )
    
    axis(side=1, at=c(0,4), labels=NA, lwd.ticks=0)
    mtext(side=1, padj=0, at = 1:3, text = c("none", "none", "ASO 1"), line=-0.2, cex=0.3)
    mtext(side=1, padj=0, at = c(1,2.5), text = c("43\ndpi", "73\ndpi"), line=1.5, cex=0.4)
    barwidth = 0.4
    axis(side=1, line=1, at=1 + c(-1,1)*barwidth, tck=0.03, labels=NA)
    axis(side=1, line=1, at=c(2,3) + c(-1,1)*barwidth, tck=0.03, labels=NA)
  
    axis(2, las = 2)

    
    if (tis == "brain") {
      mtext(side=2, line=3, text='prion titer\nlog10 seeds/g',cex = 0.8)
      
    }else {
    }
    
    stat <- titers %>%
      group_by(group,color) %>%
      summarise(
        n=n(),
        mean = mean(titer, na.rm = TRUE),
        se = sd(titer, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>% ungroup()%>%
      mutate(
        l95 = mean - 1.96 * se,
        u95 = mean + 1.96 * se,
        x = as.numeric(factor(group))
      )
    
    
    
    suppressWarnings(with(
      subset(stat, !is.na(l95) & !is.na(u95)),
      {
        arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
               angle = 90, code = 3, length = 0.04, lwd = 1.5, col = color)
        segments(x - 0.28, mean, x + 0.28, mean, lwd = 3, col = color)
      }
    ))
    
    beeswarm(titers$titer ~ titers$group,
             pch = 21, cex = 1.2,
             pwcol =titers$color,
             add = TRUE)
    
    dpi_colors <- tapply(df$color, df$group, `[`, 1)
    
    
    mtext(paste0(tis,'\nRT-QuIC'), side = 3, line = 0.5, cex = 0.6)
    
  }
}



## SZ nfl ####

nfl_sz = ellas %>%
  filter(experiment == 'SZ',
         ella_plate == 92,
         is.na(exclude),
         gn_rs == 3) %>%
  rename(animal = animal_id, expt = experiment) %>%
  left_join(blind %>% select(expt, animal, cohort, treatment), by = c('expt', 'animal')) %>%
  mutate(
    mean_conc = as.numeric(mean_conc),
    group = case_when(
      sample_dpi == 43 & treatment == 'none'  ~ '43dpi none',
      sample_dpi == 73 & treatment == 'none'  ~ '73dpi none',
      sample_dpi == 73 & treatment == 'ASO 1' ~ '73dpi ASO 1',
    ),
    color = case_when(
      group == '43dpi none'  ~ '#0001CD',
      group == '73dpi none'  ~ '#fb8500',
      group == '73dpi ASO 1' ~ '#CFB53B',
    )
  ) %>%
  filter(!is.na(group))


ref_nfl_sz <- nfl_sz$mean_conc[nfl_sz$group == '73dpi none']
nfl_sz_ks <- nfl_sz %>%
  group_by(group) %>%
  summarize(
    p_value = tryCatch(
      ks.test(mean_conc, ref_nfl_sz)$p.value,
      error = function(e) NA_real_
    ),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )


nfl_sz_stat <- nfl_sz %>%
  group_by(group, color) %>%
  summarise(
    n    = n(),
    mean = mean(mean_conc, na.rm = TRUE),
    se   = sd(mean_conc, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  mutate(
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
    x   = as.numeric(factor(group, levels = c('43dpi none', '73dpi none', '73dpi ASO 1')))
  ) %>%
  inner_join(nfl_sz_ks, by='group')



write_supp_table(nfl_sz %>% select(-color), 'SZ NfL individual')
write_supp_table(nfl_sz_stat %>% select(-color, -x), 'SZ NfL summarized')

SZ_nfl <- function() {
  par(mar = c(4, 4, 2, 1))

  ylims <- c(0, 2700)

  nfl_sz$group <- factor(nfl_sz$group, levels = c('43dpi none', '73dpi none', '73dpi ASO 1'))

  plot(
    NA, NA,
    xlim = c(0.5, 3.5),
    ylim = ylims,
    xlab = '', ylab = '',
    xaxt = 'n', yaxt = 'n',
    xaxs = 'i', yaxs = 'i',
    bty  = 'n'
  )

  
  axis(side=1, at=c(0,4), labels=NA, lwd.ticks=0)
  mtext(side=1, padj=0, at = 1:3, text = c("none", "none", "ASO 1"), line=-0.2, cex=0.3)
  mtext(side=1, padj=0, at = c(1,2.5), text = c("43\ndpi", "73\ndpi"), line=1.5, cex=0.4)
  barwidth = 0.4
  axis(side=1, line=1, at=1 + c(-1,1)*barwidth, tck=0.03, labels=NA)
  axis(side=1, line=1, at=c(2,3) + c(-1,1)*barwidth, tck=0.03, labels=NA)
  
  axis(2, las = 2)

  mtext(side = 2, line = 3, text = 'NfL (pg/ml)', cex = 0.8)

  with(
    subset(nfl_sz_stat, !is.na(l95) & !is.na(u95)),
    {
      arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
             angle = 90, code = 3, length = 0.04, lwd = 1.5, col = color)
      segments(x - 0.28, mean, x + 0.28, mean, lwd = 3, col = color)
    }
  )

  beeswarm(nfl_sz$mean_conc ~ nfl_sz$group,
           pch = 21, cex = 1.2,
           pwcol = nfl_sz$color,
           add = TRUE)
  
  

  mtext('plasma\nNfL', side = 3, line = 0.5, cex = 0.6)
}



### hamster weight varience ####
data <- read_tsv("data/weight_variance_analysis/data/raw_data/combined_data_processed.tsv", show_col_types = FALSE)

data <- data %>%
  filter((species == "mouse" & dpi < 120) | (species == "hamster" & dpi < 60)
         | (species == "rat" & dpi < 250 & dpi > 70)) %>%
  mutate(inoculated = ifelse(species == "rat", FALSE, inoculated)) %>%
  arrange(animal, dpi)

martingale_df <- data %>%
  group_by(animal) %>%
  filter(n() >= 3) %>%
  group_modify(function(df, key) {
    if(nrow(df) >= 6) {
      model <- lm(weight ~ poly(dpi, 2), data = df)
    } else {
      model <- lm(weight ~ dpi, data = df)
    }
    df %>% mutate(
      expected_weight = predict(model, newdata = df),
      martingale_val = weight / expected_weight
    )
  }) %>%
  ungroup()

martingale_increments <- martingale_df %>%
  group_by(animal) %>%
  mutate(
    delta_M = martingale_val - lag(martingale_val),
    dt = dpi - lag(dpi)
  ) %>%
  filter(!is.na(delta_M)) %>% 
  ungroup()

volatility_params <- martingale_increments %>%
  filter(inoculated == FALSE) %>% 
  group_by(species) %>%
  summarise(
    n_animals = n_distinct(animal),
    var_increments = var(delta_M, na.rm=TRUE),
    avg_dt = mean(dt, na.rm=TRUE),
    sigma_sq_daily = var_increments / avg_dt,
    sigma_daily = sqrt(sigma_sq_daily),
    .groups = "drop"
  ) %>% ungroup()


epsilon <- 0.15 # 15% weight loss

risk_table_detailed <- volatility_params %>%
  crossing(day_point = c(25, 50, 75, 100)) %>%
  mutate(
    # Kolmogorov Bound (Upper Risk Limit)
    risk_upper = pmin(1, 0.5 * (day_point * sigma_sq_daily) / (epsilon^2)),
    
    # Gaussian Estimate (Lower Risk Estimate)
    current_volatility = sigma_daily * sqrt(day_point),
    risk_lower = 2 * pnorm(-epsilon / current_volatility) # <--- FIXED: Assigned to 'risk_lower'
  ) %>%
  select(species, sigma_daily, day_point, risk_lower, risk_upper) %>%
  pivot_wider(
    names_from = day_point,
    values_from = c(risk_lower, risk_upper),
    names_glue = "day_{day_point}_{.value}"
  )

display_table_clean <- risk_table_detailed %>%
  mutate(species = str_to_title(species)) %>%
  mutate(across(starts_with("day"), ~ .x * 100)) %>% 
  transmute(
    Species = species,
    `Daily CV` = sprintf("%.1f%%", sigma_daily * 100),
    `Day 25` = sprintf("%.1f%% (\u2264 %.1f%%)", day_25_risk_lower, day_25_risk_upper)
  )

write_supp_table(display_table_clean, 'weight variance summary')





#  TgVole ####


## DATA PREP ####

make_sf_mri = function(survdata,meta=meta_mri) {
  sf = survfit(Surv(dpi, acm) ~ cohort_disp, data=survdata)
  sf$cohort = gsub('cohort_disp=','',names(sf$strata))
  sf$color = meta$color[match(sf$cohort, meta$cohort_disp)]
  return(sf)
}

master_mri=read_tsv('data/analytic/master_mri.tsv', col_types=cols()) %>% mutate(animal = as.character(animal))
weight_mri=read_tsv('data/analytic/weights_mri.tsv', col_types=cols())  %>% mutate(animal = as.character(animal))
surgery_mri=read_tsv('data/analytic/surgery_mri.tsv', col_types=cols())  %>% mutate(animal = as.character(animal))


blind_mri= read_tsv('data/analytic/blind_mri.tsv', col_types=cols())%>%
  mutate(
    dose_ug = sapply(dose_ug, function(x) {
      if (is.null(x) || length(x) == 0) return(NA_real_)
      as.numeric(x[1])
    })
  )%>%
  rename(txdpi=takedown_timepoint) %>% mutate(animal = as.character(animal))


meta_mri <- read_tsv('data/analytic/meta_mri.tsv', col_types=cols())



###  SURV ####
survdata_mri<- master_mri %>% 
  left_join(select(blind_mri,expt, animal, cohort,treatment, txdpi), by=c('animal','expt')) %>%
  mutate(death_date = as.Date(death_date),
         prion_endpoint=as.logical(prion_endpoint)) %>%
  select(expt,sex,animal,cohort, txdpi, treatment, dob, inoculation_date,death_date, prion_endpoint, acm,endpoint_comments,acm_comments) %>%
  left_join(meta_mri %>% distinct(expt,color,cohort,cohort_disp), by=c('expt','cohort'))


survdata_mri$dpi <- as.integer(survdata_mri$death_date- survdata_mri$dob)
# HU2 is always expressed in dpi, not age
is_hu2 <- survdata_mri$expt == 'HU2'
survdata_mri$dpi[is_hu2] <- as.integer(survdata_mri$death_date[is_hu2] - survdata_mri$inoculation_date[is_hu2])

###weight ####


weight_plot_mri <- weight_mri %>%
  mutate(date = as.Date(date)) %>%
  inner_join(master_mri %>% 
               select(expt, animal, dob,inoculation_date, prion_endpoint, acm),by=c('animal','expt') )%>%
  mutate(inoculation_date = as.Date(inoculation_date)) %>%
  inner_join(blind_mri %>% select(expt,cohort,animal,txdpi), by=c('animal','expt'))%>%
  inner_join(meta_mri %>% distinct(expt,color,cohort,cohort_disp), by=c('expt','cohort'))%>%
  filter(!str_detect(cohort, "-died|excluded"))%>%
  mutate(wtdpi = as.integer(date - if_else(expt == 'HU2', inoculation_date, dob))) %>%
  group_by(expt, animal) %>%
  mutate(
    baseline = first(weight),
    rel = (weight - baseline) / baseline
  ) %>%
  ungroup()%>%
  distinct(animal,date,weight,expt,dob,inoculation_date,prion_endpoint,acm,cohort,txdpi,color,cohort_disp,wtdpi,baseline,rel)



## VONH ####
### survival #### 
# vs. PrP KO


surv_ex <- survdata_mri %>%
  filter(expt %in% c('VONH')) 


write_supp_table(surv_ex %>% select(-color,-cohort), 'VONH cohort individual survival table')



surv_smry<- surv_table(surv_ex, ref_group = "PrP KO")

write_supp_table(surv_smry, 'VONH survival summary')


VONH_survival<- function(){

  par(mar=c(4,4,2,2))

  surv_ex <- survdata_mri %>%
    filter(expt %in% c('VONH'))

  sf = make_sf_mri(surv_ex,meta=meta_mri%>%filter(expt %in% c('VONH')) )

  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)), xlab='days of age')

}



### weights #### 



dat_ex<- weight_plot_mri %>%
  filter(expt %in% c('VONH')) %>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  

write_supp_table(dat_ex%>%select(-color,-cohort), 'VONH weights individual')



rel_weights_smry<- dat_ex %>%
  mutate(approx_dpi = round(wtdpi / 7) * 7)%>%
  group_by(approx_dpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>%
  ungroup() %>%
  filter(n > 1) %>%   select(approx_dpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'VONH weights summarized')



weight_summary <- make_weight_summary(weight_plot_mri%>%
                                        filter(expt == 'VONH'), ref_group = "PrP KO")

write_supp_table(weight_summary, 'VONH weights peak')

VONH_weight<-function(){
  
  
  par(mar=c(4,4,2,2))
  
  
  dat_ex<- weight_plot_mri %>%
    filter(expt %in% c('VONH')) %>%
    select(animal,cohort, wtdpi, color, rel,cohort_disp)    
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(110, 260),ylims = c(-0.25,0.4),leg=F)
  
  legend <- unique(dat_ex[c("cohort_disp","color")])
  legend <- legend[order(legend$cohort_disp), ]  
  
  legend(x='topright', y=1, bty='n', cex=.8, lwd=2, legend=legend$cohort_disp,
         col=legend$color,
         text.col=legend$color)
  
  
  
  
}


### nfl #### 

nfl_vonh= read_tsv('data/analytic/nfl.tsv', col_types=cols())%>%
  filter(experiment=='VONH') %>% 
  clean_names() %>%
  mutate(animal = as.character(animal_id)) %>%
  select(-animal_id)


plot_nfl <- nfl_vonh %>%
  select(-genotype,-treatment)%>%
  rename(expt=experiment,day=age_d)%>%
  inner_join(master_mri,by=c('animal','expt'))%>%
  inner_join(blind_mri,by=c('animal','expt','sex'))%>%
  inner_join(meta_mri,by=c('cohort','expt'))%>%
  select(cohort_disp, animal,sex,timepoint=day, nfl=plasma_nf_l_pg_ml,color) 

comb_meta = tibble(cohort_disp = c('RML','uninoculated'),
                   color = c( '#FF2017', '#5f6289'))

n_nfl = read_tsv('data/analytic/N_nfl.tsv', col_types=cols()) %>%
  select(cohort_disp=cohort, animal, timepoint, nfl=nfl_pgml) %>%
  filter(!is.na(nfl))%>%
  mutate(sex="")

n_nfl<- n_nfl%>%
  left_join(comb_meta, by='cohort_disp')

comb_nfl = rbind(n_nfl, plot_nfl) %>%
  mutate(timepoint_rounded = 30*round(timepoint/30,0))%>%
  filter(cohort_disp !='PrP KO')%>%
  arrange(cohort_disp)

write_supp_table(comb_nfl%>%select(-color), 'VONH NfL individual')

plot_nfl_vonh_smry <- comb_nfl %>%
  group_by(timepoint,cohort_disp,color) %>%
  summarise(
    n     = n(),
    mean_nfl = mean(nfl, na.rm = TRUE),
    se_nfl  = sd(nfl, na.rm = TRUE) / sqrt(n),
    l95   = mean_nfl - 1.96 * se_nfl,
    u95   = mean_nfl + 1.96 * se_nfl  ,
    .groups='keep') %>%
  ungroup()

write_supp_table(plot_nfl_vonh_smry%>%select(-color), 'VONH NfL summarized')




VONH_NFL<- function(){
  
  par(mar=c(4,2,2,3))
  
  comb_nfl <- comb_nfl %>% filter(cohort_disp!='PrP KO') %>% filter(timepoint_rounded <= 150)
  
  xlims <-c(0,160)
  ylims =c(20,23000)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',log="y")
  
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,30), lwd=0, labels=seq(0,600,30), line=-0.5)
  axis(side=1, at=seq(0,600,30), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='TgVole age (days) / RML dpi',cex=0.8)
  yats=rep(1:9,4) * rep(10^(1:4),each=9)
  ybigs=10^(1:4)
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=2, line=3, text='NfL (pg/ml)',cex=0.8)
  
  par(xpd=T)
  points(comb_nfl$timepoint_rounded,comb_nfl$nfl,col=comb_nfl$color)
  par(xpd=F)
  
  ks_vonh <- comb_nfl %>%
    group_by(timepoint_rounded, cohort_disp) %>%
    summarize(
      p_value = tryCatch(
        ks.test(nfl, comb_nfl$nfl[comb_nfl$cohort_disp == 'uninoculated' &
                                    comb_nfl$timepoint_rounded == first(timepoint_rounded)])$p.value,
        error = function(e) NA_real_
      ),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(
      sig = case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )

  plot_nfl_vonh_smry <- comb_nfl %>%
    group_by(timepoint_rounded,cohort_disp,color) %>%
    summarise(
      .groups='keep',
      n     = n(),
      mean_nfl = mean(nfl, na.rm = TRUE),
      se_nfl  = sd(nfl, na.rm = TRUE) / sqrt(n),
      l95   = pmax(min(ylims), lower(nfl)),
      u95   = pmin(max(ylims),upper(nfl))
      ) %>%
    filter(cohort_disp!='PrP KO') %>%
    filter(n > 1) %>%
    filter(timepoint_rounded <= 150) %>%
    ungroup() %>%
    left_join(ks_vonh %>% select(timepoint_rounded, cohort_disp, p_value, sig),
              by = c('timepoint_rounded', 'cohort_disp'))
  

  
  
  for (co in unique(plot_nfl_vonh_smry$cohort_disp)) {
    subs <- plot_nfl_vonh_smry %>%
      filter(cohort_disp==co)%>%
      filter(n>1)%>%
      filter(!is.na(l95) & !is.na(u95)) %>%
      arrange(timepoint_rounded)
    
    polygon(x=c(subs$timepoint_rounded,rev(subs$timepoint_rounded)),y=c(subs$u95,rev(subs$l95)),col=alpha(subs$color[1],0.5),border=NA) 
    lines(subs$timepoint_rounded, subs$mean_nfl, lwd = 2.5, col = subs$color)
    points(subs$timepoint_rounded, subs$mean_nfl, pch = 16, col = subs$color)
    
  }
  
  
  
  par(xpd=T)
  
  
  legend(90,50000,
         legend = c( "Tg(Vole)", "WT RML", "WT uninoculated") ,
         col    = c( "#f18701","#FF2017","#5f6289"),
         text.col = c( "#f18701","#FF2017","#5f6289"),
         pch    = 20,
         lwd    = 1,
         cex    = 0.8,
         lty    = 1,
         bty    = "n")
  par(xpd=F)
  

}


VONH_NFL_sex <- function(){

  sex_colors <- c('F' = '#FFC0CB', 'M' = '#87CEEB')

  par(mar=c(4,2,2,3))

  vonh_sex <- comb_nfl %>%
    filter(sex %in% c('F', 'M')) %>%
    filter(timepoint_rounded <= 150) %>%
    mutate(color = sex_colors[sex])

  xlims <- c(0, 160)
  ylims <- c(20, 23000)

  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', log='y')

  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,30), lwd=0, labels=seq(0,600,30), line=-0.5)
  axis(side=1, at=seq(0,600,30), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='TgVole age (days)', cex=0.8)

  yats <- rep(1:9,4) * rep(10^(1:4), each=9)
  ybigs <- 10^(1:4)

  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)

  mtext(side=2, line=3, text='NfL (pg/ml)', cex=0.8)

  par(xpd=T)
  points(vonh_sex$timepoint_rounded, vonh_sex$nfl, col=vonh_sex$color)
  par(xpd=F)

  vonh_sex_smry <- vonh_sex %>%
    group_by(timepoint_rounded, sex, color) %>%
    summarise(
      .groups = 'keep',
      n        = n(),
      mean_nfl = mean(nfl, na.rm=TRUE),
      se_nfl   = sd(nfl, na.rm=TRUE) / sqrt(n),
      l95      = pmax(min(ylims), lower(nfl)),
      u95      = pmin(max(ylims), upper(nfl))
    ) %>%
    filter(n > 1) %>%
    ungroup() %>%
    arrange(timepoint_rounded)

  for (s in c('F', 'M')) {
    subs <- vonh_sex_smry %>%
      filter(sex == s) %>%
      filter(!is.na(l95) & !is.na(u95))

    polygon(x=c(subs$timepoint_rounded, rev(subs$timepoint_rounded)),
            y=c(subs$u95, rev(subs$l95)),
            col=alpha(sex_colors[s], 0.5), border=NA)
    lines(subs$timepoint_rounded, subs$mean_nfl, lwd=2.5, col=sex_colors[s])
    points(subs$timepoint_rounded, subs$mean_nfl, pch=16, col=sex_colors[s])
  }

  par(xpd=T)
  legend(x='topleft', bty='n', cex=0.8, lwd=2,
         legend=c('Female', 'Male'),
         col=sex_colors[c('F', 'M')],
         text.col=sex_colors[c('F', 'M')])
  par(xpd=F)

  m = lm(nfl ~ sex + timepoint, data = comb_nfl %>% filter(cohort_disp=='Tg(Vole)'))
  summary(m)$coefficients %>%
    as_tibble(rownames = 'variable') %>%
    clean_names() %>%
    rename(p = pr_t) -> vonh_nfl_sex_model_coefficients
  write_supp_table(vonh_nfl_sex_model_coefficients, 'Linear model of TgVole NfL by sex.')

}




## VOSS ####
### RT-QuIC ####



vortq = read_tsv('data/analytic/dyr_rtquics.tsv', col_types=cols())%>%
  filter(experiment=='VOSS',
         is.na(exclude))%>%
  clean_names()%>% 
  mutate(animal_id=unlist(animal_id),
         sample=unlist(sample))%>%
  rename(animal=animal_id,expt=experiment)%>%
  inner_join(master_mri %>% select(animal,dob,death_date,expt), by=c('animal','expt'), relationship='many-to-many')%>%
  inner_join(select(blind_mri, cohort,animal,expt), by=c('animal','expt'))%>%
  inner_join(select(meta_mri,cohort_disp,color,cohort,expt), by=c('cohort','expt'))%>%
  mutate(days=as.integer(death_date-dob))



vortq_sk <- vortq %>%
  group_by(sample, days, dilution,cohort_disp,color) %>%
  summarize(
    positive_replicates = sum(positive_replicates),
    total_replicates    = sum(total_replicates),
    .groups = "drop"
  ) %>% ungroup() %>%
  arrange(sample, days, desc(dilution)) %>%
  group_by(sample, days,cohort_disp,color) %>%
  summarize(
    sk = spearman_karber(
      x = positive_replicates,
      n = total_replicates,
      d = log10(dilution)
    ),
    .groups = "drop"
  ) %>% ungroup() %>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2], llq=sk[,3] + -log10(bh_unit_factor), ulq=sk[,4] + -log10(bh_unit_factor)) %>%
  select(-sk) %>%
  ungroup() 


vortq_sk -> vortq_plot


write_supp_table(vortq_plot%>%select(-color), 'VOSS RT-quic summarized')


vortq_smry <- vortq_plot %>%select(-color) %>%
  mutate( 
    days= case_when(
      days %in% c(30, 31) ~ 30,
      days %in% c(148,150) ~ 150,
      days %in% c(171,180) ~ 180,
      TRUE ~ days
    ))%>%
  group_by(days) %>%
  summarise(
    n     = n(),
    mean_titer = mean(titer, na.rm = TRUE),
    se_titer   = sd(titer, na.rm = TRUE) / sqrt(n),
    l95   = mean_titer - 1.96 * se_titer,
    u95   = mean_titer + 1.96 * se_titer,
    .groups = "drop"
  )%>% ungroup()

write_supp_table(vortq_smry, 'VOSS RT-quic summarized by days')


VOSS_RT_quic<- function(){
  
  par(mar=c(4,4,2,3))
  
  
  
  
  par(xpd=F)
  xlims <-c(0,max(vortq_plot$days))
  ylims =c(5,11)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',log="y")
  
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,60), lwd=0, labels=seq(0,600,60), line=-0.5)
  axis(side=1, at=seq(0,600,60), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='age (days)',cex=0.8)
  yats=log10(rep(1:9,7) * rep(10^(5:11),each=9))
  ybigs=5:11
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=2, line=2, text='titer (log10 seeds/g)',cex=0.8)
  
  par(xpd=T)
  points(vortq_plot$days,vortq_plot$titer,col=vortq_plot$color)
  par(xpd=F)
  
  loess_model = loess(titer ~ days, data=vortq_plot, span=2)
  x = seq(0:180)
  y = predict(loess_model, x)
  points(x, y, type='l', col=vortq_plot$color[1], lwd=2)
  abline(h=vortq_plot$llq[1], lty=3, col='red')
  mtext(side=4, at=vortq_plot$llq[1], col='red', las=2, text='LLQ',cex=0.8)
  
  legend(x='topleft', bty='n', cex=0.8, lwd=2, legend=unique(vortq_plot$cohort_disp),
         col=unique(vortq_plot$color),
         text.col=unique(vortq_plot$color))
  
  
}


##VOAG ####
### survival ####

surv_ex <- survdata_mri %>%
  filter(expt %in% c('VOAG')) %>%
  filter(!is.na(sex),
         cohort_disp!='PrP KO' )


surv_ex$censored = is.na(surv_ex$death_date)
surv_ex$acm[surv_ex$censored] = FALSE
surv_ex$acm[!surv_ex$censored] = TRUE


write_supp_table(surv_ex %>% select(-color,-cohort), 'VOAG survival individual')



sf = survfit(Surv(dpi, acm) ~ sex, data=surv_ex)
sf$cohort = gsub('sex=','',names(sf$strata))
sf$color =c('pink','skyblue')




# survfit medians
km_med <- as.data.frame(summary(sf)$table) %>%
  tibble::rownames_to_column("sex") %>%
  mutate(
    sex = sub("sex=", "", sex),
    
  ) %>%
  select(sex, median)


# N, events, mean±sd[acm] 
basic <- surv_ex %>%
  group_by(sex) %>%
  summarise(
    N = n(),
    experiments = paste(unique(expt), collapse = ", "),
    events = sum(acm == 1),
    mean_acm = mean(dpi[acm == 1], na.rm = TRUE),
    sd_acm = sd(dpi[acm == 1], na.rm = TRUE),
    .groups = "drop"
  )%>% ungroup()

# reference group median
ref_median <- km_med$median[km_med$sex == 'F']

# percent delta from reference
km_med <- km_med %>%
  mutate(delta = (median - ref_median) / ref_median)

all_groups <- sort(unique(surv_ex$sex))
compare_groups <- setdiff(all_groups, 'F')
ref_group='F'
pairwise <- map_dfr(compare_groups, function(g) {
  
  dat2 <- surv_ex %>% 
    filter(sex %in% c(ref_group, g)) %>%
    mutate(tmp = factor(sex, levels = c(ref_group, g)))
  
  s <- suppressWarnings(survdiff(Surv(dpi, acm) ~ tmp, data = dat2))
  p <- 1 - pchisq(s$chisq, df = 1)

  tibble(sex = g, p_vs_ref = p)
}) %>%
  bind_rows(tibble(sex = ref_group, p_vs_ref = 1))


final <- basic %>%
  left_join(km_med, by = "sex") %>%
  left_join(pairwise, by = "sex") %>%
  mutate(
    median = round(median, 1),
    delta = round(delta, 2),
    mean_acm = round(mean_acm, 1),
    sd_acm = round(sd_acm, 2),
    p_vs_ref = signif(p_vs_ref, 3)
  )%>%   mutate( median = ifelse(
    is.na(median),
    paste0(">", max(surv_ex$dpi, na.rm = TRUE)), # if no median reched then use last dpi
    median
  ))%>%
  mutate(across(everything(), ~ ifelse(is.na(.x) | is.nan(.x), "", .x)))

write_supp_table(final, 'VOAG survival summary')

stat_registry$surv[[length(stat_registry$surv) + 1]] <- surv_stats(surv_ex %>% mutate(cohort_disp = sex), 'F') %>%
  mutate(experiment = 'VOAG', species = 'mouse', genotype = 'Tg(Vole)', prion_strain = 'none',
         treatment = 'none', dose_ug = '', treatment_start_dpi = '', time_axis = 'age (days)')





VOAG_survival<- function(){
  
  par(mar=c(4,4,2,0))
  surv_ex <- survdata_mri %>%
    filter(expt %in% c('VOAG')) %>%
    filter(!is.na(sex),
           cohort_disp!='KO' )
  
  
  surv_ex$censored = is.na(surv_ex$death_date)
  surv_ex$acm[surv_ex$censored] = FALSE
  surv_ex$acm[!surv_ex$censored] = TRUE
  
  
  
  sf = survfit(Surv(dpi, acm) ~ sex, data=surv_ex)
  sf$cohort = gsub('sex=','',names(sf$strata))
  
  sf$color=c('pink','skyblue')
  
  
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi,na.rm = T)),leg=F,xlab='age (days)')
  
  legend(x='topright', y=1, bty='n', cex=1, lwd=2, legend=c('Female','Male'),
         col=c('pink','skyblue'),
         text.col=c('pink','skyblue'))
  
  
  
}






## VOTE ####
### elisa ####

elisa_VOTE =  elisas %>%
  filter(is.na(exclude),
         experiment=='VOTE')%>%
  mutate(animal_id=unlist(animal_id),
         sample=unlist(sample),
         ngml_av=unlist(ngml_av),
         se_mean=unlist(se_mean))%>%
  rename(animal=animal_id,expt=experiment)

MoPosHiQC <- elisas %>%
  filter(experiment == "VOTE",
         sample == "Mo Pos Hi QC") %>%
  mutate(ngml_av = unlist(ngml_av)) %>%
  pull(ngml_av)




plot_elisa_VOTE <- elisa_VOTE %>%
  left_join(blind_mri %>% select(animal, sex,cohort,genotype, treatment, inoculum,expt), by=c('animal','expt'))%>%
  left_join(meta_mri%>% select(cohort,cohort_disp,color),by='cohort')%>%
  mutate(
    color_sex = case_when(
      sex %in% c("F") ~ "pink",   
      sex %in% c("M")   ~ "lightblue"
    ),
    moposhiqc = MoPosHiQC
  ) %>%
  mutate(aso_mean = mean(ngml_av[cohort_disp=='None'], na.rm=T)) %>%
  mutate(norm_conc = ngml_av / aso_mean,
         fold = ngml_av / moposhiqc )

write_supp_table(plot_elisa_VOTE%>%select(-color_sex,-color,-exclude,-cohort), 'VOTE ELISA individual')


t_male_female <- plot_elisa_VOTE %>%
  summarize(
    p_value = tryCatch(
      t.test(norm_conc ~ sex)$p.value,
      error = function(e) NA_real_
    )
  )%>%
  mutate(test = "female_vs_male")




t_ASO_vs_PBS <- plot_elisa_VOTE %>%
  group_by(cohort_disp) %>%
  summarize(
    p_value = if (unique(cohort_disp) == "None") {
      NA_real_
    } else {
      tryCatch(
        t.test(
          plot_elisa_VOTE$norm_conc[plot_elisa_VOTE$cohort_disp == "None"],
          plot_elisa_VOTE$norm_conc[plot_elisa_VOTE$cohort_disp == unique(cohort_disp)]
        )$p.value,
        error = function(e) NA_real_
      )
    }
  )%>% ungroup()%>%filter(!is.na(p_value))%>%
  mutate(test = paste0(cohort_disp, "_vs_None"))


combined_VOTE_t_test <- t_male_female %>%
  select(test, p_value) %>%
  bind_rows(t_ASO_vs_PBS %>% select(test, p_value))

write_supp_table(combined_VOTE_t_test, 'VOTE ELISA T-test')




LLQ <- 10 / aso_mean * 100

stat_save <- plot_elisa_VOTE %>%
  group_by(cohort_disp,color) %>%
  summarise(
    mean = mean(norm_conc, na.rm = TRUE),
    mean_fold = mean(fold, na.rm = TRUE),
    se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )%>% ungroup() %>%
  mutate(
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
    x = as.numeric(factor(cohort_disp, levels = c("None", "ASO 6 300 µg","ASO 6 500 µg")))
  )

write_supp_table(stat_save%>%select(-color,-x), 'VOTE ELISA summarized')


stat_save_2 <- plot_elisa_VOTE %>%
  group_by(cohort_disp,sex) %>%
  summarise(
    mean = mean(norm_conc, na.rm = TRUE),
    se = sd(norm_conc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>% ungroup()%>%
  mutate(
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
  )



write_supp_table(stat_save_2, 'VOTE ELISA summarized by sex and treatment')



VOTE_elisa<- function(){
  par(mar=c(4,5,2,5))
  
  par(xpd=F)
  n_by_group <- plot_elisa_VOTE %>% count(cohort_disp, name = "n")
  x_label <- setNames(
    paste0(levels(factor(plot_elisa_VOTE$cohort_disp, levels = c("None", "ASO 6 300 µg","ASO 6 500 µg"))),
           "\n(n=", n_by_group$n[match(c("None", "ASO 6 300 µg","ASO 6 500 µg"), n_by_group$cohort_disp)], ")"),
    c("None", "ASO 6 300 µg","ASO 6 500 µg")
  )
  
  plot_elisa_VOTE$cohort_disp <- factor(plot_elisa_VOTE$cohort_disp, levels = c("None", "ASO 6 300 µg","ASO 6 500 µg"))
  
  plot(
    NA, NA,
    xlim = c(0.5, 3.5),
    ylim = c(0, 1.6),
    xlab = "", ylab = "residual PrP",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i",
    bty  = "n"
  )
  
  axis(1, at = 1:3, labels = FALSE, tick = FALSE)
  
  axis(1, at = 1:3, labels = x_label[c("None", "ASO 6 300 µg","ASO 6 500 µg")], line=0,col = NA,  tick = FALSE ,cex.axis=0.8)
  axis(1, at=0:4, labels=NA, tck=FALSE)
  
  axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2) # in terms of norm_conc
  right_axis_factor = plot_elisa_VOTE$aso_mean[1]/plot_elisa_VOTE$moposhiqc[1]
  axis(side=4, at=0:10/right_axis_factor, labels=0:10, las=2) # in terms of fold
  
  par(xpd=T)
  
  text(4,0.9,
       "Fold WT",
       srt = -90,
       adj = c(-0.2, 0.5),
       cex = 0.9
  )
  par(xpd=F)
  
  
  
  stat <- plot_elisa_VOTE %>%
    group_by(cohort_disp,color) %>%
    summarise(
      mean = mean(norm_conc, na.rm = TRUE),
      l95 = lower(norm_conc),
      u95 = upper(norm_conc),
      mean_fold = mean(fold, na.rm = TRUE),
      l95_fold = lower(fold),
      u95_fold = upper(fold),
      sexpval = t.test(norm_conc[sex=='F'],norm_conc[sex=='M'])$p.value,
      sexratio= mean(norm_conc[sex=='F'])/mean(norm_conc[sex=='M']),
      .groups = "drop"
    ) %>% ungroup()%>%
    mutate(
      x = as.numeric(factor(cohort_disp, levels = c("None", "ASO 6 300 µg","ASO 6 500 µg")))
    )
  write_supp_table(stat, 'TgVole ASO target engagement data summarized.')
  
  with(stat, {
    arrows(x0 = x, y0 = l95, x1 = x, y1 = u95,
           angle = 90, code = 3, length = 0.08, lwd = 1.5,col ='grey')
    segments(x - 0.25, mean, x + 0.25, mean, lwd = 3,col = 'grey')
  })
  
  
  plot_elisa_VOTE <- plot_elisa_VOTE %>%
    arrange(desc(norm_conc))
  
  
  beeswarm(plot_elisa_VOTE$norm_conc ~ plot_elisa_VOTE$cohort_disp,
           pch =21, cex = 1.2,
           pwcol =plot_elisa_VOTE$color_sex,
           add = TRUE)
  
  
  cols <- plot_elisa_VOTE %>%
    distinct(cohort_disp, color,color_sex,sex) %>%
    arrange(cohort_disp)
  
  
}













## VOE1 ####
### survival####

surv_ex <- survdata_mri %>%
  filter(expt %in% c("VOE1"))


write_supp_table(surv_ex %>% select(-color,-cohort), 'VOE1 survival individual')



surv_smry<- surv_table(surv_ex, ref_group = "PBS")%>%
  arrange(match(cohort_disp, c('PBS','ASO 6 300 µg')))


write_supp_table(surv_smry, 'VOE1 survival summarized')


VOE1_survival<- function(){
  par(mar=c(4,4,2,1))
  
  surv_ex <- survdata_mri %>%
    filter(expt %in% c("VOE1"))
  
  
  
  anno_date <- surgery_mri %>% filter(expt %in% c('VOE1'))%>%filter(event_id=='icv')%>%
    inner_join(surv_ex,by=c('animal','expt'))%>%
    mutate(anno_dpi=as.integer(date-dob))%>%
    pull(anno_dpi)%>%
    unique(
    )
  
  anno <- c(60,150)
  
  
  
  
  sf = make_sf_mri(surv_ex)
  
  
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
  
  points(
    x = anno,
    y = rep(1.04, 2),
    pch = 25,                     
    bg = '#30336b',        
    col = '#30336b',      
    cex = 0.5
  )
  
  
}





### weight ####

dat_ex<- weight_plot_mri %>%
  filter(expt == 'VOE1') %>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  

write_supp_table(dat_ex%>% select(-color), 'VOE1 weights individual')

rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            se = sd(rel, na.rm = TRUE) / sqrt(n()),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup() %>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'VOE1 weights summarized')


weight_summary <- make_weight_summary(weight_plot_mri%>%
                                        filter(expt == 'VOE1'), ref_group = "PBS")%>%
  arrange(match(cohort_disp, c('PBS','ASO 6 300ug')))

write_supp_table(weight_summary, 'VOE1 weights peak')



VOE1_weight<-function(){
  
  par(mar=c(4,5,2,4))
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(50, 160),ylims = c(-0.2,0.4))
  
  mtext('VOE1', side = 3, line = 0.5, font = 2, cex = 0.9)
}
### nfl####

nfl_VOE1 = read_tsv('data/analytic/nfl.tsv', col_types=cols())%>%
  filter(experiment=='VOE1') %>%
  clean_names() %>%
  mutate(animal = as.character(animal_id))



plot_nfl_VOE1<- nfl_VOE1 %>%
  mutate(treatment=unlist(treatment))%>%
  select(-genotype,-treatment)%>%
  rename(expt=experiment,day=age_d)%>%
  left_join(blind_mri %>% select(expt,cohort,animal,genotype), by=c('expt','animal'))%>%
  left_join(meta_mri %>% distinct(expt,cohort,color,cohort_disp), by=c('expt','cohort'))%>%
  filter(!is.na(cohort))%>%
  mutate(cohort_disp=factor(cohort_disp,levels=c('PBS','ASO 6 300 µg')))%>%
  rename(nfl=plasma_nf_l_pg_ml)

write_supp_table(plot_nfl_VOE1%>%select(-color,-cohort), 'VOE1 NfL individual')


plot_nfl_voe1_smry <- plot_nfl_VOE1%>%select(-cohort) %>%
  group_by(day,cohort_disp,color) %>%
  summarise(
    n     = n(),
    mean_nfl = mean(nfl, na.rm = TRUE),
    se_nfl  = sd(nfl, na.rm = TRUE) / sqrt(n),
    l95   = mean_nfl - 1.96 * se_nfl,
    u95   = mean_nfl + 1.96 * se_nfl  ,
    .groups='drop')%>% ungroup()

write_supp_table(plot_nfl_voe1_smry%>%select(-color), 'VOE1 NfL summarized')



VOE1_NFL<- function(){
  
  par(mar=c(4,4,2,2))
  
  
  
  
  
  par(xpd=F)
  xlims <-c(0,max(plot_nfl_VOE1$day))
  ylims =c(20,20000)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',log="y")
  
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,60), lwd=0, labels=seq(0,600,60), line=-0.5)
  axis(side=1, at=seq(0,600,60), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='Age(days)',cex=0.8)
  yats=rep(1:9,4) * rep(10^(1:4),each=9)
  ybigs=10^(1:4)
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=2, line=3, text='NfL (pg/ml)',cex=0.8)
  
  par(xpd=T)
  points(plot_nfl_VOE1$day,plot_nfl_VOE1$nfl,col=plot_nfl_VOE1$color)
  par(xpd=F)
  
  for (co in unique(plot_nfl_voe1_smry$cohort_disp)) {
    subs <- plot_nfl_voe1_smry %>%
      filter(cohort_disp==co)%>%
      filter(n>1)
    
    polygon(x=c(subs$day,rev(subs$day)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,0.5),border=NA)
    lines(subs$day, subs$mean_nfl, lwd = 2.5, col = subs$color)
    points(subs$day, subs$mean_nfl, pch = 16, col = subs$color)
    
  }
  
  legend(x='bottomright', y=1, bty='n', cex=.8, lwd=2, legend=unique(plot_nfl_voe1_smry$cohort_disp),
         col=unique(plot_nfl_voe1_smry$color),
         text.col=unique(plot_nfl_voe1_smry$color))
  
  
  
}




## VOE2 ####
### survival####

surv_ex <- survdata_mri %>%
  filter(expt %in% c("VOE2"))%>%
  mutate(cohort_disp=factor(cohort_disp,levels=c('PrP KO','PBS','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg')))%>%
  arrange(cohort_disp)

write_supp_table(surv_ex%>%select(-cohort,-color), 'VOE2 survival table')



surv_smry<- surv_table(surv_ex, ref_group = "PBS")%>%
  arrange(match(cohort_disp, c('PBS','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg','PrP KO')))



write_supp_table(surv_smry, 'VOE2 survival summary')

VOE2_survival<- function(){
  
  surv_ex <- survdata_mri %>%
    filter(expt %in% c("VOE2"))%>%
    mutate(cohort_disp=factor(cohort_disp,levels=c('PBS','PrP KO','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg')))%>%
    arrange(cohort_disp)
  
  anno_date <- surgery_mri %>% filter(expt %in% c('VOE2'))%>%
    inner_join(surv_ex,by=c('animal','expt'))%>%
    mutate(anno_dpi=as.integer(date-dob))%>%
    pull(anno_dpi)%>%
    unique(
    )
  
  inj1=round(mean(anno_date[seq(1, length(anno_date), by = 2)]), 0)
  inj2=round(mean(anno_date[seq(2, length(anno_date), by = 2)]), 0)
  
  
  par(mar=c(4,2,2,2))
  sf = make_sf_mri(surv_ex,meta = meta_mri%>% filter(expt=="VOE2"))
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
  
  points(
    x = c(inj1,inj2),
    y = rep(1.04, 2),
    pch = 25,                     
    bg = 'black',        
    col = 'black',      
    cex = 0.5
  )
  
}






### weight ####

dat_ex<- weight_plot_mri %>%
  filter(expt == 'VOE2') %>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)  %>%
  mutate(cohort_disp=factor(cohort_disp,levels=c('PBS','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg','PrP KO')))%>%
  arrange(cohort_disp)



Wednesday <-  as.Date("2023-02-08") + seq(0, 624, 7)

Wednesday_mapping <- tibble(
  that_Wednesday = rep(Wednesday, each = 7),
  date        = rep(Wednesday, each = 7) + rep(0:6, times = length(Wednesday))
)

voe2_indiv <-  weight_plot_mri %>%
  filter(expt == 'VOE2')%>%
  inner_join(Wednesday_mapping, by = "date") %>%
  mutate(
    Wednesday_dpi       = as.numeric(as.Date(that_Wednesday) - as.Date(dob)),
    Wednesday_wpi       = Wednesday_dpi %/% 7,
    Wednesday_approx_dpi = Wednesday_wpi * 7
  ) %>%
  select(animal,  weight,approx_dpi=Wednesday_approx_dpi,rel,cohort_disp,color,wtdpi)%>%
  arrange(animal,approx_dpi)


write_supp_table(voe2_indiv%>% select(-color), 'VOE2 weights individual')

rel_weights_smry<- voe2_indiv %>%
  group_by(approx_dpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup() %>%
  filter(n > 1) %>%   select(approx_dpi,cohort_disp,n,mean,l95,u95,color)



rel_weights_smry_save<- weight_plot_mri %>%
  filter(expt == 'VOE2')%>%
  inner_join(Wednesday_mapping, by = "date") %>%
  mutate(
    Wednesday_dpi       = as.numeric(as.Date(that_Wednesday) - as.Date(dob)),
    Wednesday_wpi       = Wednesday_dpi %/% 7,
    Wednesday_approx_dpi = Wednesday_wpi * 7,
    approx_dpi=Wednesday_approx_dpi
  ) %>%
  group_by(approx_dpi,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n()))%>% ungroup() %>%
  filter(n > 1) %>%   select(approx_dpi,cohort_disp,n,mean,l95,u95)


write_supp_table(rel_weights_smry_save, 'VOE2 weights summarized')


weight_summary <- make_weight_summary(voe2_indiv%>% mutate(expt = 'VOE2'), ref_group = "PBS")%>%
  arrange(match(cohort_disp, c('PBS','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg','PrP KO')))



write_supp_table(weight_summary, 'VOE2 weights peak')


VOE2_weight<-function(){
  
  par(xpd=F)
  par(mar=c(4,4,2,2))
  dat_ex<- voe2_indiv %>%
    select(animal, wtdpi=approx_dpi, color, rel,cohort_disp)   %>%
    mutate(cohort_disp=factor(cohort_disp,levels=c('PBS','PrP KO','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg')))%>%
    arrange(cohort_disp)%>%
    filter(cohort_disp != 'PrP KO') 
  
  
  
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(0, 370),ylims = c(-0.25,0.4),leg = F)
  
  
  legend <- unique(dat_ex[c("cohort_disp","color")])
  legend <- legend[order(legend$cohort_disp), ]  
  
  legend(x=100, y=0.56, xpd=NA, ncol=2, bty='n', cex=.7, lwd=2, legend=legend$cohort_disp,
         col=legend$color,
         text.col=legend$color)
  
  
  
  
}


### NfL####


nfl_VOE2 = read_tsv('data/analytic/nfl.tsv', col_types=cols())%>%
  filter(experiment=='VOE2') %>%
  clean_names() %>% 
  mutate(animal = as.character(animal_id)) %>%
  select(-animal_id)





plot_nfl_VOE2<- nfl_VOE2 %>%
  mutate(treatment=unlist(treatment))%>%
  rename(expt=experiment,day=age_d)%>%
  rename(nfl=plasma_nf_l_pg_ml)%>%
  inner_join(blind_mri %>% select(animal,cohort,expt),by=c('animal','expt'))%>%
  inner_join(meta_mri %>% select(expt,cohort,cohort_disp,color),by=c('cohort','expt'))%>%
  mutate(cohort_disp=factor(cohort_disp,levels=c('PBS','PrP KO','ASO 6 100 µg','ASO 6 300 µg','ASO 6 500 µg')))


write_supp_table(plot_nfl_VOE2%>% select(-cohort,-color), 'VOE2 NfL individual')



plot_nfl_voe2_smry <- plot_nfl_VOE2%>%select(-cohort) %>%
  group_by(day,cohort_disp,color) %>%
  summarise(
    n     = n(),
    mean_nfl = mean(nfl, na.rm = TRUE),
    se_nfl  = sd(nfl, na.rm = TRUE) / sqrt(n),
    l95   = mean_nfl - 1.96 * se_nfl,
    u95   = mean_nfl + 1.96 * se_nfl  ,
    .groups='drop')%>% ungroup()

write_supp_table(plot_nfl_voe2_smry%>%select(-color), 'VOE2 NfL summarized')




VOE2_NFL<- function(){
  
  
  plot_nfl_VOE2<- plot_nfl_VOE2%>%
    filter(cohort_disp != 'PrP KO') 
  
  par(mar=c(4,4,2,2))
  par(xpd=F)
  xlims <-c(0,max(plot_nfl_VOE2$day))
  ylims =c(20,20000)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',log="y")
  
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,60), lwd=0, labels=seq(0,600,60), line=-0.5)
  axis(side=1, at=seq(0,600,60), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='Age(days)',cex=0.8)
  yats=rep(1:9,4) * rep(10^(1:4),each=9)
  ybigs=10^(1:4)
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=2, line=2.75, text='NfL (pg/ml)',cex=0.8)
  
  par(xpd=T)
  points(plot_nfl_VOE2$day,plot_nfl_VOE2$nfl,col=plot_nfl_VOE2$color)
  par(xpd=F)
  
  
  plot_nfl_voe2_smry<- plot_nfl_voe2_smry%>%
    filter(cohort_disp != 'PrP KO') 
  
  for (co in unique(plot_nfl_voe2_smry$cohort_disp)) {
    subs <- plot_nfl_voe2_smry %>%
      filter(cohort_disp==co)%>%
      filter(n>1)%>%
      mutate( l95 = pmax(l95, ylims[1]))
    
    
    polygon(x=c(subs$day,rev(subs$day)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,0.5),border=NA)
    lines(subs$day, subs$mean_nfl, lwd = 2.5, col = subs$color)
    points(subs$day, subs$mean_nfl, pch = 16, col = subs$color)
    
  }
  
  legend(x='bottomright', y=1, bty='n', cex=.8, lwd=2, legend=unique(plot_nfl_voe2_smry$cohort_disp),
         col=unique(plot_nfl_voe2_smry$color),
         text.col=unique(plot_nfl_voe2_smry$color))
  
  
}




# Humanized ####
## HA  ####

### survival ####


surv_ex <- survdata %>%
  filter(expt %in% c("HA"))

write_supp_table(surv_ex%>%select(-color), 'HA survival individual')


sf = make_sf(surv_ex)

surv_smry<- surv_table(surv_ex, ref_group = "uninoculated")%>%
  arrange(match(cohort_disp, c('uninoculated','E200K MM','E200K MV','MM1 sCJD isolate 1',
                               'MM1 sCJD isolate 2','V2 sCJD isolate 1','VV2 sCJD isolate 2')))


write_supp_table(surv_smry, 'HA survival summarized')

HA_HB_legend<- function(){
  
  surv_ex <- survdata %>%
    filter(expt %in% c("HA"))%>%
    mutate(
      cohort_disp = factor(
        cohort_disp,
        levels = c("uninoculated","MM1 sCJD isolate 1", "MM1 sCJD isolate 2", "VV2 sCJD isolate 1", "VV2 sCJD isolate 2","E200K MM","E200K MV")
      )
    )
  
  sf = make_sf(surv_ex)
  
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  
  legend(x='center', y=1, bty='n', cex=.6, lwd=2, legend= sf$cohort,
         col=sf$color,
         text.col=sf$color)
  
  
}

HA_survival<- function(){
  
  
  par(mar=c(4,4,2,2))
  
  surv_ex <- survdata %>%
    filter(expt %in% c("HA"))%>%
    mutate(
      cohort_disp = factor(
        cohort_disp,
        levels = c("uninoculated","MM1 sCJD isolate 1", "MM1 sCJD isolate 2", "VV2 sCJD isolate 1", "VV2 sCJD isolate 2","E200K MM","E200K MV")
      )
    )
  
  sf = make_sf(surv_ex)
  
  
  
  xlims <- c(0, max(surv_ex$dpi))
  ylims = c(0,1.05)
  plot(sf, col=sf$color, lwd=2, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
  mtext(side=2, line=2.75, text='survival',cex=0.8)
  par(xpd=T)
  
  
  
}

## HB  ####
### survival ####


surv_ex <- survdata %>%
  filter(expt %in% c("HB"),cohort_disp !='exclude')

write_supp_table(surv_ex %>% select(-color) , 'HB survival individuals')


surv_smry<- surv_table(surv_ex, ref_group = "uninoculated")%>%
  arrange(match(cohort_disp, c('uninoculated','E200K MM','E200K MV','MM1 sCJD isolate 1',
                               'MM1 sCJD isolate 2','V2 sCJD isolate 1','VV2 sCJD isolate 2')))


write_supp_table(surv_smry, 'HB survival summarized')

HB_survival<- function(){
  par(mar=c(4,2,2,2))
  
  surv_ex <- survdata %>%
    filter(expt %in% c("HB"),cohort_disp !='exclude')%>%
    mutate(
      cohort_disp = factor(
        cohort_disp,
        levels = c("uninoculated","MM1 sCJD isolate 1", "MM1 sCJD isolate 2", "VV2 sCJD isolate 1", "VV2 sCJD isolate 2","E200K MM","E200K MV")
      )
    )
  
  sf = make_sf(surv_ex)
  
  
  par(xpd=F)
  xlims <- c(0, max(surv_ex$dpi))
  ylims = c(0,1.05)
  plot(sf, col=sf$color, lwd=2, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
  mtext(side=2, line=2.75, text='survival',cex=0.8)
  
  
  
  
} 

### HA/HB RT-QuIC ####

rtq_save <- bind_rows(
  `00061` = read.table("data/analytic/00061_summary.tsv", sep = "\t", header = TRUE),
  `00062` = read.table("data/analytic/00062_summary.tsv", sep = "\t", header = TRUE),
  `00063` = read.table("data/analytic/00063_summary.tsv", sep = "\t", header = TRUE),
  .id = "source"
)%>%
  select(animal=sample, everything())%>%
  left_join(blind %>% select(expt,animal,genotype,cohort,inoculum,treatment,dose_ug,treatment_start_dpi),by=c('animal'))%>%
  left_join(master %>% select(animal,dob,death_date),by=c('animal'))%>%
  left_join(meta %>% select(cohort,cohort_disp,color),by='cohort')%>%
  select(source,expt,animal,cohort_disp,genotype,dob,death_date,inoculum,treatment,dose_ug,treatment_start_dpi,color,everything())



write_supp_table(rtq_save%>%filter(expt=='HA')%>%select(-cohort,-inoculum,-color), 'HA RT-quic individual')

rtq_save %>%
  left_join(master %>%
              select(-cage) %>%
              mutate(dpi = as.integer(death_date - inoculation_date)), by=c('animal','expt'))%>%
  group_by(expt, animal, cohort_disp, genotype, dpi, treatment,prion_endpoint, endpoint_comments) %>%
  summarize(.groups='keep',
            sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution)))%>% ungroup()%>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor),
         se=sk[,2]) %>%
  select(-sk) %>% ungroup() -> titers_ha_hk

write_supp_table(titers_ha_hk%>%filter(!is.na(expt),expt=='HA'), 'HA RT-quic summarized')











rtq <- bind_rows(
  `00061_avkinetic` = read.table("data/analytic/00061_avkinetic.tsv", sep = "\t", header = TRUE),
  `00062_avkinetic` = read.table("data/analytic/00062_avkinetic.tsv", sep = "\t", header = TRUE),
  .id = "source"
)






plot_rt_quic_ha_hb <-  rtq_save %>%
  select(expt, animal, genotype, color, treatment, inoculum=cohort_disp, dilution, positive_replicates, total_replicates) %>%
  filter(!is.na(expt))


plot_rt_quic_ha_hb %>%
  mutate(
    genotype= case_when(
      expt == "HA" ~ "Tg26372 hom",
      expt == "HK" ~ "Ki817 hom",
      TRUE ~ genotype))%>%
  group_by(expt, animal, genotype, inoculum,treatment, color) %>%
  summarize(.groups='keep',
            sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution))) %>% ungroup()%>%
  mutate(titer = sk[,1] + -log10(bh_unit_factor), se=sk[,2]) %>%
  select(-sk) %>% ungroup() -> titers_ha_hb


titers_ha_hb <- titers_ha_hb %>% arrange(desc(titer))%>%
  mutate(
    label = paste0(inoculum, " → ", genotype),
    l95 = titer - 1.96 * se,
    u95 = titer + 1.96 * se
  )

HA_HB_HK_rt_quic_all <- function(){ 
  
  
  df <- titers_ha_hb%>%
    mutate(
      label = factor(
        label,
        levels = c("E200K MM → Tg26372 hom","E200K MV → Tg26372 hom", "MM1 sCJD isolate 1 → Tg26372 hom","MM1 sCJD isolate 2 → Tg26372 hom", "MM1 sCJD isolate 2 → Ki817 hom", "VV2 sCJD isolate 1 → Ki817 hom")
      )
    )%>%arrange(desc(label))
  
  par(mar=c(4,12,2,2))  
  par(xpd=F)
  plot(
    NA, NA,
    xlim = c(6,12),
    ylim = c(0, 8.4),
    axes = FALSE,
    ann = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  mtext(side = 1, line = 2, text = "titer (log10 seeds/g)",cex =0.8)
  axis(side = 1, at = 0:13, labels = NA, tck = -0.025)
  axis(side = 1, at = 0:13, tck = -0.025)
  
  rect(
    xleft   = 0,
    ybottom = seq_len(nrow(df)) - 0.4,
    xright  = df$titer,
    ytop    = seq_len(nrow(df)) + 0.4,
    col     = adjustcolor(df$color, alpha.f = 0.3),
    border  = NA
  )
  
  
  par(xpd=T)
  
  suppressWarnings(arrows(df$l95, seq_len(nrow(df)),
         df$u95, seq_len(nrow(df)),
         angle = 90, length = 0.05,
         code = 3,       
         col = df$color, lwd = 2))
  par(xpd=F)
  
  segments(df$l95, seq_len(nrow(df)),
           df$u95, seq_len(nrow(df)),
           col=df$color, lwd=2)
  
  axis(2, at=0:9 ,labels=NA, tck=FALSE)
  
  axis(2,
       at = seq_len(nrow(df)),
       labels = df$label,
       las = 2, tick = FALSE, cex.axis = 0.8)

  # Flag the anomalous bar: HAG3 (E200K MV → Tg26372 hom) died of an
  # intercurrent cause before reaching terminal prion disease, so its titer
  # is unrepresentative of the cohort. Mark its bar with an asterisk.
  anom <- which(df$animal == "HAG3")
  if (length(anom) > 0) {
    par(xpd = TRUE)
    text(x = df$u95[anom] + 0.25, y = anom,
         labels = "*", cex = 1.6, font = 2)
    par(xpd = FALSE)
  }

}









HA_HB_HK_rt_quic_indiv<- function(s){
  
  par(mar=c(4,4,3,3))
  
  
  
  
  sample_id <- sub("_1e.*$", "", rownames(rtq))
  
  
  
  hrs <- suppressWarnings(as.numeric(gsub("X", "", colnames(rtq))))
  maxhr_index <- max(which(hrs <= 23.5))
  
  
  idx <- which(sample_id == s)
  
  plot(
    NA, NA,
    xlim = c(0, 24),
    ylim = c(-0.05, 1.05),
    axes = FALSE,
    ann = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  
  axis(side = 1, at = 0:24, labels = NA, tck = -0.025)
  axis(side = 1, at = 0:4 * 6, tck = -0.05)
  mtext(side = 1, line = 2, text = "time (h)",cex =0.8)
  
  axis(side = 2, at = 0:4 / 4, labels = NA, las = 2)
  
  if (s %in% c("HAA1", "HKF1")) {
    
    axis(side = 2, at = 0:4 / 4, labels = percent(0:4 / 4), las = 2)
    
    mtext(side = 2, line = 3, text = "relative fluorescence", cex = 0.8)
  }
  
  
  
  for (i in idx) {
    
    dilution <- as.numeric(sub("^.*_", "", rownames(rtq)[i]))
    alpha_val <- (12 + log10(dilution)) / 12
    fadecolor <- alpha("#a9a9a9", alpha_val)

    points(hrs[!is.na(hrs)], rtq[i, !is.na(hrs)], type = "l", lwd = 1, col = fadecolor)
    
    mtext(
      side = 4,
      at = rtq[i, maxhr_index],
      las = 2,
      line = 0.25,
      col = fadecolor,
      text = substr(rownames(rtq)[i], 6, 11),
      cex=0.5
    )
    
    
    
    i = which(grepl('RES-04',rownames(rtq)))
    dilution = substr(rownames(rtq)[i],10,11)
    fadecolor = alpha('#eb4d4b', (12-as.numeric(dilution))/12)
    points(hrs[!is.na(hrs)], rtq[i, !is.na(hrs)], type='l', lwd=1, col=fadecolor)
    
    if (s %in% c("HAA1", "HKF1")) {
      mtext(side=4, at=rtq[i,maxhr_index], las=2, line=0.25, col=fadecolor, text='poscon',cex=0.5)
    }
    i = which(grepl('RES-11',rownames(rtq)))
    dilution = substr(rownames(rtq)[i],10,11)
    fadecolor = alpha('#3a6ea5', (12-as.numeric(dilution))/12)
    points(hrs[!is.na(hrs)], rtq[i, !is.na(hrs)], type='l', lwd=1, col=fadecolor)
    
    if (s %in% c("HAA1", "HKF1")) {
      mtext(side=4, at=rtq[i,maxhr_index], las=2, line=0.25, col=fadecolor, text='negcon',,cex=0.5)
      
    }
    
  }
  
  
  
}


## HG/HH  ####
### survival ####


surv_ex <- survdata %>%
  filter(expt %in% c("HG","HH"))%>%
  mutate(
    surv_group = case_when(
      str_detect(cohort_disp, regex("none",    ignore_case = TRUE)) ~ "None",
      TRUE ~ cohort_disp
    ))


write_supp_table(surv_ex %>% select(-color,-cohort,-surv_group) , 'HG/HH survival individual')


surv_smry <- surv_table(surv_ex, ref_group = "Tg26372 hom (20 copies)")%>%
  arrange(match(cohort_disp, c('Tg26372 hom (20 copies)','uninoculated','PrP KO','Tg25109 het (3 copies)','Tg26372 het (10 copies)')))


write_supp_table(surv_smry, 'HG/HH survival individual')


HG_HH_survival<- function(){
  
  par(mar=c(4,4,2,2))
  surv_ex <- survdata %>%
    filter(expt %in% c("HG","HH"))
  
  sf=make_sf(surv_ex)
  
  
  
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)))
  
  
  
  
  
  
  
}

### weights ####


dat_ex<- weight_plot %>%
  filter(expt %in% c("HG","HH"))%>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp)


write_supp_table(dat_ex%>%select(-color,-cohort), 'HG/HH weights individual')

rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'HG/HH weights summarized')





weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt %in% c("HG","HH")), ref_group = 'Tg26372 hom (20 copies)')%>%
  arrange(match(cohort_disp, c('Tg26372 hom (20 copies)','uninoculated','PrP KO','Tg25109 het (3 copies)','Tg26372 het (10 copies)')))


write_supp_table(weight_summary, 'HG/HH weights peak')





HG_HH_weight<- function(){
  
  par(mar=c(4,3,2,2))
  
  
  
  dat_ex<- weight_plot %>%
    filter(expt %in% c("HG","HH"))%>%
    select(cohort, weight,wtdpi, color, rel,cohort_disp)
  
  
  
  
  
  rel_weights_smry<- dat_ex %>%
    group_by(wtdpi, color,cohort_disp) %>%
    summarize(.groups='keep',
              mean = mean(rel),
              sd = sd(rel),
              n = n(),
              l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
              u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
    filter(n > 1) 
  
  xlims = c(120, 400)
  ylims = c(-0.25,0.24)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=seq(0,1000,50), labels=NA, tck=-0.05)
  axis(side=1, at=seq(0,1000,100), lwd=0, labels=seq(0,1000,100), line=-0.5)
  axis(side=1, at=seq(0,1000,10), labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation')
  axis(side=2, at=seq(-2,15,.2), labels=NA, las=2, tck=-0.05)
  axis(side=2, at=seq(-2,15,.2), labels=percent(seq(-2,15,.2), signed=T), lwd=0, las=2, line=-0.5)
  axis(side=2, at=seq(-2,15,.1), labels=NA, las=2, tck=-0.025)
  abline(h=0)
  abline(h=-.2, col='red', lty=3)
  text(
    x = xlims[2]-10,     
    y = -0.15,
    labels = "endpoint",
    col = "red",
    adj = c(1, 1),
    cex = 0.8
  ) 
  mtext(side=2, line=2.75, text='weight (% change)',cex=0.8)
  
  for (this_cohort in unique(rel_weights_smry$cohort_disp)) {
    subs = rel_weights_smry %>% filter(cohort_disp==this_cohort) 
    points(subs$wtdpi, subs$mean, type='l', lwd=2, col=subs$color)
    polygon(x=c(subs$wtdpi,rev(subs$wtdpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,ci_alpha),border=NA)
  }
  
  
  
  
}

### nests ####


dat_ex <- nest_plot %>%
  filter(expt %in% c("HG","HH"))%>%
  filter(!is.na(cohort))

write_supp_table(dat_ex%>%select(-color,-cohort), 'HG/HH nests individual')


HG_HH_nests<- function(){
  par(mar=c(4,5,2,4))
  
  
  plot_nest(dat_ex,
            xlims = c(0, 500))
  
  mtext('HG/HH', side = 3, line = 0.5, font = 2, cex = 0.9)
  
}
## HK  ####
### survival ####

surv_ex <- survdata %>%
  filter(expt %in% c("HK"))



write_supp_table(surv_ex %>% select(-color) , 'HK survival table')


surv_smry<- surv_table(surv_ex, ref_group = "uninoculated")%>%
  arrange(match(cohort_disp, c('uninoculated','E200K MM','E200K MV','MM1 sCJD isolate 1','MM1 sCJD isolate 2','VV2 sCJD isolate 1','VV2 sCJD isolate 2')))


write_supp_table(surv_smry, 'HK survival summary')


HK_survival<- function(){
  
  par(mar=c(4,4,2,2))
  
  surv_ex <- survdata %>%
    filter(expt %in% c("HK"))%>% 
    mutate(
      cohort_disp = factor(
        cohort_disp,
        levels = c("uninoculated","MM1 sCJD isolate 1", "MM1 sCJD isolate 2", "VV2 sCJD isolate 1", "VV2 sCJD isolate 2","E200K MM","E200K MV")
      )
    )
  
  
  
  sf = make_sf(surv_ex)
  plot_sf(sf,xlims = c(0, max(surv_ex$dpi)),leg=F)
}



### RT-QuIC ####

rtq_63 = read.table('data/analytic/00063_avkinetic.tsv',sep='\t',header=T)%>%
  mutate(source='00063_avkinetic')






write_supp_table(rtq_save%>%filter(expt=='HK')%>%select(-cohort), 'HK RT- quic individual')

write_supp_table(titers_ha_hk%>%filter(!is.na(expt),expt=='HK'), 'HK RT-quic summarized')














## HJ ####
### survival ####

surv_ex <- survdata %>%
  filter(expt %in% c("HJ"))


write_supp_table(surv_ex %>% select(-color,-cohort) , 'HJ survival table')


surv_smry<- surv_table(surv_ex, ref_group = "Ki817 / Ki817 VV2 sCJD")%>%
  arrange(match(cohort_disp, c('Ki817 / Ki817 VV2 sCJD','Ki817 / ZH3 VV2 sCJD','uninoculated')))


write_supp_table(surv_smry, 'HJ survival summary')


HJ_survival<- function(){
  
  par(mar=c(4,4,2,1))
  
  surv_ex <- survdata %>%
    filter(expt %in% c("HJ"))
  sf = make_sf(surv_ex)
  
  plot_sf(sf,xlims = c(0, 600))
  
  
}





### weights ####

dat_ex<- weight_plot %>%
  filter(expt %in% c("HJ"))%>%
  select(animal,cohort, weight,wtdpi, color, rel,cohort_disp) 

write_supp_table(dat_ex%>%select(-color,-cohort), 'HJ weights individual')

rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'HJ weights summarized')




weight_summary <- make_weight_summary(weight_plot%>%
                                        filter(expt %in% c("HJ")), ref_group = "Ki817 / Ki817 VV2 sCJD")%>%
  arrange(match(cohort_disp, c('Ki817 / Ki817 VV2 sCJD','Ki817 / KO VV2 sCJD','uninoculated')))


write_supp_table(weight_summary, 'HG/HH weights peak')





HJ_weight<- function(){
  
  
  par(mar=c(4,4,2,2))
  dat_ex<- weight_plot %>%
    filter(expt %in% c("HJ"))%>%
    select(animal,cohort, weight,wtdpi, color, rel,cohort_disp) 
  
  
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(40, 600),ylims = c(-0.25,0.6),leg=FALSE)
  
  meta %>%
    filter(expt=='HJ') %>%
    distinct(cohort_disp, color) -> leg
  par(xpd=T)
  legend(x=130,y=0.75,leg$cohort_disp,col=leg$color,lwd=2,bty='n',cex=0.8)
  par(xpd=F)
  
}
### nests ####


dat_ex <- nest_plot %>%
  filter(expt %in% c("HJ"))%>%
  filter(!is.na(cohort)) 


write_supp_table(dat_ex%>%select(-color,-cohort), 'HJ nests individual')

HJ_nest<- function(){
  
  par(mar=c(4,5,2,4))
  
  
  plot_nest(dat_ex,
            xlims = c(0, 500))
  
  mtext('HJ', side = 3, line = 0.5, font = 2, cex = 0.9)
  
}


## HU2 ####
### survival ####

surv_ex <- survdata_mri %>%
  filter(expt %in% c("HU2"))


write_supp_table(surv_ex %>% select(-color) , 'HU2 survival table')


surv_smry<- surv_table(surv_ex, ref_group = "uninoculated")%>%
  arrange(match(cohort_disp, c('uninoculated','RES-05','RES-09')))

write_supp_table(surv_smry, 'HU2 survival summary')


HU2_survival<- function(){
  
  par(mar=c(4,4,2,1))
  surv_ex <- survdata_mri %>%
    filter(expt %in% c("HU2"))
  
  
  
  sf = make_sf_mri(surv_ex,meta = meta_mri%>% filter(expt=="HU2"))
  
  
  plot_sf(sf,xlims = c(0, 676))
}


### weights ####


dat_ex<- weight_plot_mri %>%
  filter(expt == 'HU2') %>%
  select(animal,weight, cohort, wtdpi, color, rel,cohort_disp)     

write_supp_table(dat_ex%>%select(-color,-cohort), 'HU2 weights individual')


rel_weights_smry<- dat_ex %>%
  group_by(wtdpi, color,cohort_disp) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>% ungroup()%>%
  filter(n > 1) %>%   select(wtdpi,cohort_disp,n,mean,l95,u95,color)

write_supp_table(rel_weights_smry%>%ungroup() %>%select(-color), 'HU2 weights summarized')



weight_summary <- make_weight_summary(weight_plot_mri%>%
                                        filter(expt %in% c("HU2")), ref_group = "uninoculated")%>%
  arrange(match(cohort_disp, c('uninoculated','RES-05','RES-09')))

write_supp_table(weight_summary, 'HU2 weights peak')







HU2_weight<- function(){
  
  par(mar=c(4,4,2,1))
  
  dat_ex<- weight_plot_mri %>%
    filter(expt == 'HU2') %>%
    select(animal,weight, cohort, wtdpi, color, rel,cohort_disp)     
  
  plot_wts(weights_indiv = dat_ex,
           xlims = c(50, 700),ylims = c(-0.2,0.5),leg=F)
  
}


### NfL ####




nfl_hu2 = ellas %>%
  filter(experiment=='HU2',
         is.na(exclude))%>%
  mutate(animal=as.character(animal_id),
         expt=experiment)%>%
  select(-animal_id) %>%
  left_join(blind_mri %>% select(expt,cohort,animal,genotype), by=c('expt','animal'))%>%
  left_join(meta_mri %>% distinct(expt,cohort,color,cohort_disp), by=c('expt','cohort'))%>%
  select(animal,color,species,expt,tissue,sample_dpi,genotype,cohort,cohort_disp,sample_name,mean_rfu,rfu_percent_cv,mean_conc,conc_percent_cv,cohort_disp)





plot_nfl_hu2<- nfl_hu2 %>%
  mutate(animal=as.numeric(animal))%>%
  filter(!is.na(animal))%>%
  group_by(expt, sample_dpi,color,genotype,cohort_disp,cohort) %>%
  summarise(
    mean = mean(mean_conc, na.rm = TRUE),
    sd   = mean(mean_conc, na.rm = TRUE),
    n = n(),
    l95 = mean- 1.96*sd(mean_conc)/sqrt(n()),
    u95 = mean + 1.96*sd(mean_conc)/sqrt(n()),
    .groups   = "drop"
  )%>% ungroup()

write_supp_table(plot_nfl_hu2%>%select(-color,-cohort), 'HU2 NfL summarized')


ks_by_dpi_hu2<- nfl_hu2 %>%
  filter(cohort_disp %in% c("sCJD VV2 isolate 1", "uninoculated")) %>%
  group_by(sample_dpi) %>%
  summarize(
    n_sCJD_VV2   = sum(cohort_disp == "sCJD VV2 isolate 1"),
    n_uninoc = sum(cohort_disp == "uninoculated"),
    p_value = if (n_sCJD_VV2 >= 2 & n_uninoc >= 2) {
      tryCatch(
        ks.test(
          mean_conc[cohort_disp == "sCJD VV2 isolate 1"],
          mean_conc[cohort_disp == "uninoculated"]
        )$p.value,
        error = function(e) NA
      )
    } else {
      NA
    },
    .groups = "drop"
  )%>% ungroup()%>%
  mutate(
    sig = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

write_supp_table(ks_by_dpi_hu2%>%select(-sig), 'Hu2 NfL ks test')

HU2_nfl <- function(){ 
  par(mar=c(4,4,2,2))
  par(xpd=F)
  xlims <-c(180,max(plot_nfl_hu2$sample_dpi))
  ylims =c(20,2000)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',log="y")
  
  axis(side=1, at=seq(0,600,10), labels=NA, tck=-0.025)
  axis(side=1, at=seq(0,600,60), lwd=0, labels=seq(0,600,60), line=-0.5)
  axis(side=1, at=seq(0,600,60), labels=NA, tck=-0.05)
  mtext(side=1, line=1.5, text='days post-inoculation',cex=0.8)
  yats=rep(1:9,4) * rep(10^(1:4),each=9)
  ybigs=10^(1:4)
  
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=0, labels=ybigs)
  
  mtext(side=2, line=2.75, text='NfL (pg/ml)',cex=0.8)
  
  par(xpd=T)
  points(nfl_hu2$sample_dpi,nfl_hu2$mean_conc,col=plot_nfl_hu2$color)
  par(xpd=F)
  
  for (co in unique(plot_nfl_hu2$cohort)) {
    subs <- plot_nfl_hu2 %>%
      filter(cohort==co)
    
    polygon(x=c(subs$sample_dpi,rev(subs$sample_dpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,0.5),border=NA)
    lines(subs$sample_dpi, subs$mean, lwd = 2.5, col = subs$color)
    points(subs$sample_dpi, subs$mean, pch = 16, col = subs$color)
    
  }
  
  
  star_pos <- plot_nfl_hu2 %>%
    group_by(sample_dpi) %>%
    summarize(
      y = max(u95, mean, na.rm = TRUE),
      .groups = "drop"
    )%>% ungroup()
  
  
  anno <- left_join(ks_by_dpi_hu2, star_pos, by = "sample_dpi")
  
  
  text(
    x = anno$sample_dpi,
    y = anno$y +500,  
    labels = anno$sig,
    cex = 1,
    font = 2
  )
  
  
}











### PrP expression by ELISA ####

mouse_genotypes <- read_tsv('data/analytic/mouse_genotypes.tsv', col_types=cols(id=col_character()))

mouse_genotype_params <- read_tsv('data/analytic/mouse_genotypes_params.tsv', col_types=cols()) %>%
  bind_rows(tibble(
    gt    = c('TacKI/TacKI', 'TacKI/-'),
    color = c('#005B96',      '#6497B1')
  ))

qc_samples <- c('Mo Neg QC', 'Mo Pos Hi QC', 'Mo Pos Lo QC', 'Mo Pos Mid QC')

elisa_hu <- bind_rows(
  read_tsv('data/analytic/026_summary.tsv', col_types=cols()) %>% mutate(plate = 26L),
  read_tsv('data/analytic/027_summary.tsv', col_types=cols()) %>% mutate(plate = 27L)
) %>%
  filter(!sample %in% qc_samples) %>%
  rename(id = sample) %>%
  left_join(mouse_genotypes %>% select(id, gt), by = 'id') %>%
  left_join(mouse_genotype_params, by = 'gt')

wt_means_by_plate <- elisa_hu %>%
  filter(gt == 'wild-type') %>%
  group_by(plate) %>%
  summarise(wt_mean = mean(ngml_av, na.rm = TRUE), .groups = 'drop')

elisa_hu <- elisa_hu %>%
  left_join(wt_means_by_plate, by = 'plate') %>%
  mutate(rel = ngml_av / wt_mean)

plot_gt_order <- c('TacKI/-', 'TacKI/TacKI', '-/-', '+/-', 'wild-type')
plot_gt_disp  <- c('TacKI/TacKI' = 'Ki817/Ki817',
                   'TacKI/-'     = 'Ki817/ZH3',
                   '-/-'         = 'ZH3/ZH3',
                   '+/-'         = 'ZH3/+',
                   'wild-type'   = 'WT')

mashup <- elisa_hu %>%
  filter(gt %in% plot_gt_order) %>%
  mutate(y    = match(gt, plot_gt_order),
         disp = plot_gt_disp[gt])

smry <- mashup %>%
  group_by(gt, color, y, disp) %>%
  summarise(
    mean_raw = mean(ngml_av, na.rm = TRUE),
    mean     = mean(rel,     na.rm = TRUE),
    se       = sd(rel,       na.rm = TRUE) / sqrt(n()),
    .groups  = 'drop'
  ) %>%
  mutate(
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se
  )

write_supp_table(mashup %>% select(-color, -wt_mean, -y, -disp), 'HuPrP ELISA individual')
write_supp_table(smry   %>% select(-color, -y),                   'HuPrP ELISA summarized')

hu_prp_elisa <- function() {
  par(mar=c(4,8,3,1))
  dilution     <- 200
  llq          <- 0.05
  meanwt       <- smry$mean_raw[smry$gt == 'wild-type']
  llq_relative <- llq * dilution / meanwt
  barwidth     <- 0.33
  xlims <- c(0, max(c(smry$u95, mashup$rel), na.rm = TRUE) * 1.15)
  ylims <- range(smry$y) + c(-0.5, 0.5)

  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=pretty(xlims), tck=-0.05, labels=NA)
  axis(side=1, at=pretty(xlims, n=20), tck=-0.02, labels=NA)
  axis(side=1, at=pretty(xlims), line=-0.5, lwd=0)
  mtext(side=1, line=1.5, text='brain PrP (fold WT)', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=smry$y, text=smry$disp, col='#000000', las=2, line=0.25, cex=0.8)
  abline(v=llq_relative, lwd=0.75, lty=3, col='red')
  abline(v=1, lty=3)
  mtext(side=3, at=llq_relative, cex=0.75, font=2, line=0.25, text='   LLQ', col='red')
  suppressWarnings(arrows(x0=smry$l95, x1=smry$u95, y0=smry$y,
                          angle=90, length=0.05, code=3, lwd=1, col=smry$color))
  rect(xleft=rep(0, nrow(smry)), xright=smry$mean,
       ybottom=smry$y - barwidth, ytop=smry$y + barwidth,
       lwd=1, col=alpha(smry$color, ci_alpha), border=NA)
  points(x=mashup$rel, y=mashup$y, pch=21, col=mashup$color, bg='#FFFFFF')
}







# SUPPLEMENT #### 


## usinf descriptive name

expand_title <- function(title, cohort_lookup) {
  
  first_token <- sub(" .*", "", title)
  
  if (first_token %in% cohort_lookup$cohort) {
    long_name <- cohort_lookup$display[
      match(first_token, cohort_lookup$cohort)
    ]
    return(sub(first_token, long_name, title))
  }
  
  title
}

supplement_directory <- supplement_directory %>%
  mutate(title = vapply(title, expand_title, character(1), cohort_lookup = cohort_lookup))


# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)


# STANDALONE STAT TABLES (not part of the supplement workbook) ####

stat_cols <- c('experiment','species','genotype','prion_strain','treatment','dose_ug',
               'treatment_start_dpi','group','control_group','time_axis','figure_panel')

# prion isolate characterization studies are omitted from Tables S1 and S2
omit_expts <- c('HA', 'HB', 'HK')

# studies that span several experiment codes are sorted together so the reference group comes first
study_key <- function(x) case_when(x %in% c('HG', 'HH') ~ 'HG/HH',
                                   x %in% c('SD', 'SE; SK') ~ 'SC/SD/SE/SK',
                                   x %in% c('SG', 'SH') ~ 'SG/SH',
                                   TRUE ~ x)

panel_surv <- c(RR = 'Figure 1A', RG = 'Figure 1F', RP = 'Figure 1H', RE = 'Figure S2G',
                SD = 'Figure 2C', `SE; SK` = 'Figure 2C', SG = 'Figure S3B', SH = 'Figure S3B',
                VONH = 'Figure 3A', VOAG = 'Figure 3D', VOE1 = 'Figure 3F', VOE2 = 'Figure 3G',
                HG = 'Figure 4C', HH = 'Figure 4C', HJ = 'Figure 4F', HU2 = 'Figure S5E')
panel_wt   <- c(RR = 'Figure 1B', RG = 'Figure S2A', RP = 'Figure S2H',
                SD = 'Figure S3C', `SE; SK` = 'Figure S3C',
                VONH = 'Figure S4A', VOE2 = 'Figure S4D',
                HG = 'Figure 4D', HH = 'Figure 4D', HJ = 'Figure S5I', HU2 = 'Figure S5F')

table_s01 <- bind_rows(stat_registry$surv) %>%
  distinct(experiment, group, .keep_all = TRUE) %>%
  filter(!experiment %in% omit_expts) %>%
  mutate(figure_panel = unname(panel_surv[experiment])) %>%
  mutate(event_rate = paste0(events, '/', n)) %>%
  select(experiment, control_group, time_axis, n, events,
         species, genotype, prion_strain, treatment, dose_ug, treatment_start_dpi,
         group, figure_panel, event_rate, median_survival, p_vs_control) %>%
  arrange(study_key(experiment), group != control_group, group) %>%
  mutate(across(c(median_survival, p_vs_control), ~ ifelse(is.na(.x), '\u2014', as.character(.x))))
write_tsv(table_s01, 'display_items/table-s01.tsv', na = '')

table_s02 <- bind_rows(stat_registry$peak) %>%
  distinct(experiment, group, .keep_all = TRUE) %>%
  filter(!experiment %in% omit_expts) %>%
  mutate(figure_panel = unname(panel_wt[experiment])) %>%
  select(experiment, control_group, time_axis, n,
         species, genotype, prion_strain, treatment, dose_ug, treatment_start_dpi,
         group, figure_panel, mean_peak_dpi, l95, u95, p_vs_control) %>%
  mutate(across(c(mean_peak_dpi, l95, u95), ~ round(.x, 1))) %>%
  arrange(study_key(experiment), group != control_group, group) %>%
  mutate(p_vs_control = ifelse(is.na(p_vs_control), '\u2014', as.character(p_vs_control)))
write_tsv(table_s02, 'display_items/table-s02.tsv', na = '')








# FIGURES ####


### Figure 1. Rats ####
tell_user('done.\nCreating Figure 1...')
if ('figure-1'=='figure-1') {
  
  resx=600
  png("display_items/figure-1.png",width = 6.5*resx, height = 8*resx, res = resx)
  
  
  
  
  layout_matrix <- matrix(
    c(
      1, 1,1, 2, 2,      
      3, 3,3, 4, 4,     
      5, 5,6, 6, 6,     
      7, 8,9, 9,  9       
    ),
    nrow = 4,
    byrow = TRUE
  )
  
  layout(layout_matrix, widths=c(6,4,5,5,5), heights=c(4,4,4,4))
  
  
  
  RR_survival()
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  RR_weight()
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)
  
  
  RR_nfl()
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)
  RR_rtqic()
  mtext("D", side=3, cex=1.5, adj = 0, line = 0)
  
  
  rat_elisa()
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)
  RG_survival()
  mtext("F", side=3, cex=1.5, adj = 0, line = 0)
  
  
  RD_Brain_RNA()
  mtext("G", side = 3, adj = -0.6, line = 0, cex = 1.5)
  RD_CSF()
  RP_survival()
  mtext("H", side=3, cex=1.5, adj = 0, line = 0)
  
  silence_is_golden = dev.off()
  
}

### Figure S2. More about Rats ####
tell_user('done.\nCreating Figure S2...')
if ('figure-S2'=='figure-S2') {
  
  resx=600
  png("display_items/figure-s2.png",width = 6.5*resx, height = 8*resx, res = resx)
  
  
  layout_matrix <- matrix(
    c(
      1,2,
      3,4,
      5,6,
      7,7,
      8,9
      
    ),
    nrow = 5,
    byrow = TRUE
  )
  
  
  layout(layout_matrix, heights=c(1,1,1,.2,1))
  panel = 1
  
  RG_weight()
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0); panel = panel + 1

  
  regions       <- c('cortex', 'brain stem', 'hippocampus', 'thalamus', 'spinal cord')
  region_labels <- c('cortex', 'brainstem',  'hippocampus', 'thalamus', 'spinal cord')
  col_it  <- '#458B00'
  col_icv <- '#BA74FF'
  dose_centers <- setNames(c(1, 2.5, 4, 5.5), c('100', '300', '1000', '3000'))
  
  for (i in seq_along(regions)) {
    d <- rat_it_icv %>%
      filter(region == regions[i]) %>%
      mutate(
        x_center = dose_centers[as.character(dose_ug)],
        x_pos    = x_center + ifelse(route == 'IT', -0.2, 0.2),
        col      = ifelse(route == 'IT', col_it, col_icv)
      )
    
    stat <- d %>%
      group_by(dose_ug, route) %>%
      summarise(
        mean  = mean(residual, na.rm = TRUE),
        se    = sd(residual, na.rm = TRUE) / sqrt(n()),
        x_pos = first(x_pos),
        col   = first(col),
        .groups = 'drop'
      ) %>%
      mutate(
        l95 = mean - 1.96 * se,
        u95 = mean + 1.96 * se
      )
    
    xlims = c(0.5, 6.5)
    plot(NA, NA,
         xlim = xlims,
         ylim = c(0, 1.5),
         xlab = '', ylab = '',
         xaxt = 'n', yaxt = 'n',
         xaxs = 'i', yaxs = 'i',
         bty  = 'n')
    
    axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
    axis(2, at = seq(0, 1.5, by = 0.5), labels = percent(seq(0, 1.5, by = 0.5)), las = 2)
    mtext(side=2, line=3, text='residual Prnp', cex=0.8)
    axis(1, at = c(1, 2.5, 4, 5.5),
         labels = c('100', '300', '1000', '3000'),
         tick = FALSE, line = -0.5)
    mtext(expression(paste('dose (', mu, 'g)')), side = 1, line = 1.5, cex = 0.8)
    
    abline(h = 1, lty = 2, col = 'gray60')
    
    suppressWarnings(with(stat, {
      arrows(x0 = x_pos, y0 = l95, x1 = x_pos, y1 = u95,
             angle = 90, code = 3, length = 0.05, lwd = 1.5, col = col)
      segments(x_pos - 0.12, mean, x_pos + 0.12, mean, lwd = 3, col = col)
    }))
    
    set.seed(42)
    points(jitter(d$x_pos, amount = 0.05), d$residual,
           pch = 21, cex = 0.8, bg = d$col, col = d$col)
    
    mtext(region_labels[i], side = 3, line = 0.2, font = 1)
    mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0); panel = panel + 1
  }
  
  par(mar=c(0,0,0,0))
  plot.new()
  legend('center', horiz=T, legend = c('IT', 'ICV'), pch = 21,
         pt.bg = c(col_it, col_icv), col = c(col_it, col_icv),
         bty = 'n', cex = 1.2, pt.cex = 1.2)
  
  
  RE_survival()
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0); panel = panel + 1
  RP_weight()
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0); panel = panel + 1
  
  
  
  
  silence_is_golden = dev.off()
  
}


### Figure 2. Hamsters #####

tell_user('done.\nCreating Figure 2...')
if ('figure-2'=='figure-2') {
  
  
  resx=300
  png("display_items/figure-2.png",width = 6.5*resx, height = 5*resx, res = resx)
  
  
  layout_matrix <- matrix(
    c(
      1,1,2,2,3,3,3,
      4,4,4,4,6,6,6,
      5,5,5,5,7,8,9
      
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  layout(layout_matrix, widths=c(1,1,1,1,1.33,1,1.33))
  
  
  
  
  hamster_elisa()
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  hamster_aso('BH')
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)
  SC_SD_SE_SK_survival()
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)
  
  SC_SD_SE_SK_rtquic('brain')
  mtext("D", side=3, cex=1.5, adj = 0, line = 0)
  SC_SD_SE_SK_rtquic('CSF')
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)
  
  
  SC_SD_SE_SK_nfl()
  mtext("F", side=3, cex=1.5, adj = 0, line = 0)
  
  SZ_rtqic()
  mtext("G", side=3, cex=1.5, adj = -3.5, line = 0)
  mtext("H", side=3, cex=1.5, adj = -0.6, line = 0)
  
  SZ_nfl()
  mtext("I", side=3, cex=1.5, adj = 0, line = 0)
  
  silence_is_golden = dev.off()
  
}

### Figure S3. More about Hamsters ####
tell_user('done.\nCreating Figure S3...')
if ('figure-S3'=='figure-S3') {
  
  
  png("display_items/figure-s3.png",width = 6.5*resx, height = 4*resx, res = resx)
  
  layout_matrix <- matrix(
    c(
      1, 1,2, 2,
      3, 3, 3,3 
      
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  layout(layout_matrix, widths=c(5,5,5,5), heights=c(4,4))
  
  
  
  
  
  hamster_aso('CSF')
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  
  SG_SH_survival()
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)
  
  SC_SD_SE_SK_weight()
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)
  
  
  
  silence_is_golden = dev.off()
}




## Figure 3. TgVole ####


tell_user('done.\nCreating Figure 3...')
if ('figure-3'=='figure-3') {
  
  png("display_items/figure-3.png",width = 6.5*resx, height = 6*resx, res = resx)
  
  
  layout_matrix <- matrix(
    c(
      1, 1,2, 2,3,3,
      4,4,4, 5,5,5,
      6,6,6,7, 7,7
      
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  layout(layout_matrix, widths=c(5,5,5,5,5,5), heights=c(4,4,4))
  
  
  VONH_survival()
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  VONH_NFL()
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)
  VOSS_RT_quic()
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)
  
  
  VOAG_survival()
  mtext("D", side=3, cex=1.5, adj = 0, line = 0)
  VOTE_elisa()
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)
  
  
  VOE1_survival()
  mtext("F", side=3, cex=1.5, adj = 0, line = 0)
  VOE2_survival()
  mtext("G", side=3, cex=1.5, adj = 0, line = 0)
  
  
  
  silence_is_golden = dev.off()
}

### Figure S4. More about TgVole ####
tell_user('done.\nCreating Figure S4...')
if ('figure-s4'=='figure-s4') {
  
  
  
  png("display_items/figure-s4.png",width = 6.5*resx, height = 6*resx, res = resx)

  layout_matrix <- matrix(
    c(
      1, 1, 2, 2,
      3, 3, 4, 4,
      5, 5, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  layout(layout_matrix, widths=c(5,5,5,5), heights=c(4,4,4))




  VONH_weight()
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)

  VONH_NFL_sex()
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)

  VOE1_NFL()
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)

  VOE2_weight()
  mtext("D", side=3, cex=1.5, adj = 0, line = 0)

  VOE2_NFL()
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)
  
  
  silence_is_golden = dev.off()
  
}

## Figure 4. Humanized mice ####


tell_user('done.\nCreating Figure 4...')
if ('figure-4'=='figure-4') {
  
  resx=600
  png("display_items/figure-4.png",width = 6.5*resx, height = 6*resx, res = resx)
  
  
  layout_matrix <- matrix(
    c(
      7,7, 1,1,1,2,2,2,
      3,3,3,3,  4,4,4,4,
      5,5,5,5,  6,6,6,6
    ),
    nrow = 3,
    byrow = TRUE
  )

  layout(layout_matrix, widths=c(1, 1, 3.5, 2, 1, 2.5, 1.5, 2.5), heights=c(4,4,4))





  HA_survival()
  mtext(side=3, line=0.75, cex=0.7, text='Tg26372 PRNP 129M')
  mtext(side=3, line=-0.25, cex=0.6, text='Prion isolate characterization')
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  HB_survival()
  mtext(side=3, line=0.5, cex=0.7, text='Tg25109 PRNP 129M')
  mtext(side=3, line=-0.25, cex=0.6, text='Prion isolate characterization')
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)

  HG_HH_survival()
  mtext(side=3, line=0.75, cex=0.7, text='Tg PRNP 129M genetic dose-response')
  mtext(side=3, line=-0.25, cex=0.6, text='sCJD MM1 isolate 1')
  mtext("C", side=3, cex=1.5, adj = 0, line = 0)
  HG_HH_weight()
  mtext(side=3, line=0.75, cex=0.7, text='Tg PRNP 129M genetic dose-response')
  mtext(side=3, line=-0.25, cex=0.6, text='sCJD MM1 isolate 1')
  mtext("D", side=3, cex=1.5, adj = 0, line = 0)

  HK_survival()
  mtext(side=3, line=0.5, cex=0.7, text='Ki817 PRNP 129V')
  mtext(side=3, line=-0.25, cex=0.6, text='Prion isolate characterization')
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)

  HJ_survival()
  mtext(side=3, line=0.75, cex=0.7, text='Ki817 PRNP 129V genetic dose-response')
  mtext(side=3, line=-0.25, cex=0.6, text='sCJD VV2 isolate 1')
  mtext("F", side=3, cex=1.5, adj = 0, line = 0)
  
  HA_HB_legend()
  
  silence_is_golden = dev.off()
  
}


### Figure S5. More about humanized mice ####
tell_user('done.\nCreating Figure S5...')
if ('figure-S5'=='figure-S5') {
  
  
  resx=600
  png("display_items/figure-s5.png",width = 6.5*resx, height = 8.5*resx, res = resx)
  
  
  layout_matrix <- matrix(
    c(
      1,1,1, 2,2,2,2,2,2, 3,3,3,
      4,4,4,4, 4,4,4,4, 4,4,4,4, 
      5,5,5,5, 6,6,6,6, 7,7,7,7, 
      8,8,8,8, 9,9,9,9, 10,10,10,10
    ),
    nrow = 4,
    byrow = TRUE
  )

  layout(layout_matrix, widths=rep(5,12), heights=c(4,4,4,4))




  HA_HB_HK_rt_quic_indiv('HAA1')
  mtext("A", side=3, cex=1.5, adj = 0, line = 0)
  HA_HB_HK_rt_quic_all()
  mtext("B", side=3, cex=1.5, adj = 0, line = 0)
  HA_HB_HK_rt_quic_indiv('HKF1')
  mtext("C", side=3, cex=1.5, adj = -0.55, line = 0)
  
  
  par(mar=c(0.5,0,1.5,0))
  raster_panel = image_convert(image_read('data/images/ki817_vv2_histo.png'),'png')
  plot(as.raster(raster_panel))
  mtext("D", side=3, cex=1.5, adj = 0, line = -2)

  HU2_survival()
  mtext("E", side=3, cex=1.5, adj = 0, line = 0)
  HU2_weight()
  mtext("F", side=3, cex=1.5, adj = 0, line = 0)
  HU2_nfl()
  mtext("G", side=3, cex=1.5, adj = 0, line = 0)

  hu_prp_elisa()
  mtext("H", side=3, cex=1.5, adj = -0.4, line = 0)
  HJ_weight()
  mtext("I", side=3, cex=1.5, adj = 0, line = 0)
  plot.new()
  
  
  
  
  silence_is_golden = dev.off()
  
}





elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))

