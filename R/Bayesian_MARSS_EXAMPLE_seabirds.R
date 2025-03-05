library(cmdstanr)
library(tidyverse)
library(magrittr)
library(loo)
source("R/marss_stan_functions.R")

#=== Data
bird_df=read.csv("data/Seabirds.csv")

bird_df %>% group_by(Family) %>%
  summarise(n=n())

# switch to long formate for massaging
# You will have to change your dataframe's column names to match these
bird_df_long = bird_df %>% 
  # filter(Use == "yes") %>%
  # mutate(CommonName = Common.Name, Site = Location.of.population, RMU = Subspecies) %>%
  mutate(Indiv.pop.map = Site) %>% 
  ungroup() %>%
  # mutate(ts_id = paste(CommonName,"#",Site)) %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Year",
    names_prefix = "X",
    values_to = "Count",
    values_drop_na = TRUE) %>%
  mutate(Count=as.numeric(Count),Year=as.numeric(Year)) %>%
  mutate(log.spawner = log(Count+1)) %>%
  group_by(ID) %>%
  drop_na(Count) %>%
  mutate(n=n()) %>% filter(n>4) %>% # filters out short time series
  arrange(ID) %>%
  mutate(add_up = 1) # this says add up all time series in the species abundance, this is bad if you have 2 time series for 1 population
  # Below: auto sets up which time series to add up at the end, 
  # this below example is only including the largest populations of a site for an entire family, it doesnt make ecological sense but it's an example
  # group_by(Site,ID) %>%
  # mutate(med_Count = median(Count,na.rm=T)) %>%
  # group_by(Site,Family) %>%
  # mutate(add_up = case_when(med_Count==max(med_Count) ~ 1,
  #                           .default= 0))


# Necessary variables:
# ID = Time series ID, multiple IDs is possible for a single "Indiv.pop.map" or "Site"
# Indiv.pop.map = States that the time series belong to, akin to the Z matrix
# [Taxa] = Species, Family, Genus; at the taxonomic resolution you want to run the MARSS
# log.spawner = logged count data, model is in log space
# year = year of the data 
# add_up = select which time series you want to add up, generally only 1 time series per species-site has a "1", all others are "0"

# Example by filtering to Frigatebirds
species_i = "Pelecanoididae"
species_i_pre_df = bird_df_long #%>% filter(Family == species_i) 

#=== Massage data back to wide
species_i_df <- species_i_pre_df %>% 
  ungroup() %>% 
  dplyr::select(CommonName,Year,Site, log.spawner) %>%  # get just the columns that I need
  # drop_na(log.spawner) %>%
  group_by(CommonName,Site) %>% 
  complete(Year = min(Year):max(Year),
           fill = list(log.spawner = NA)) %>%
  ungroup() %>%
  pivot_wider(names_from = c("Site","CommonName"), values_from = "log.spawner") %>% 
  arrange(Year) %>% 
  column_to_rownames(var = "Year") %>% # make the years rownames
  # filter(if_any(everything(), ~ !is.na(.)))%>%
  as.matrix() %>% # turn into a matrix with year down the rows
  t() # make time across the columns

#=== Set up for Stan
# site_vec <- species_i_pre_df %>%  distinct(Site)%>% pull(Site)
# vector of populations or locations, this is your state vector
# there can be multiple time series for a single population/location
# Indiv.pop.map is the population name
id_vec <- species_i_pre_df %>% ungroup() %>% distinct(ID,Indiv.pop.map) %>% pull(Indiv.pop.map) 
id_vec <- species_i_pre_df %>% ungroup() %>% distinct(CommonName,Indiv.pop.map) %>% 
  mutate(state_vec = as.factor(row_number())) %>% pull(state_vec)

# hypotheses on states and observation time series mapping
n_states <- max(as.numeric(forcats::fct_inorder(id_vec))) # site_vec
# n_states <- max(as.numeric((id_vec))) # site_vec
r_unequal <- seq(1, nrow(species_i_df))
r_equal <- rep(1, nrow(species_i_df))
uq_unequal <- seq(1, n_states)
uq_equal <- rep(1, n_states)

#=== MARSS specifications
marss_list = list(states = as.numeric(fct_inorder(id_vec)),# as.numeric(forcats::fct_inorder(id_vec)), #site_vec
                  obsVariances = r_equal, #as.numeric(factor(id_vec)), 
                  proVariances = uq_equal,
                  trends = uq_equal,
                  stocks = uq_unequal) 

# Obs sigma by survey method?
Obs_method = F
if (Obs_method == T){
  method_vec = species_i_pre_df %>% subset(CommonName == species_i) %>% 
    distinct(ID,Sampling.method) %>%
    pull(Sampling.method)
  method_vec_in = as.numeric(fct_inorder(method_vec))
} else { method_vec_in = r_equal
method_vec = "sigma_obs"}


# set up data_list
data_list_tmp <- setup_data(y = species_i_df,
                            est_nu = FALSE,
                            est_trend = FALSE,
                            family = "gaussian", 
                            # mcmc_list = mcmc_list,
                            marss = marss_list)
data_list_tmp$data$add_up = species_i_pre_df %>% ungroup() %>% select(Site,CommonName, add_up) %>% distinct() %>% pull(add_up)
mcmc_list  = list(n_mcmc = 2500, n_burn = 1000, n_chain = 3, n_thin = 100,step_size=0.9,adapt_delta=0.8)

model_type = "Ueq"

# Let's only do Base Model
if (model_type == "Ueq") {cmd_file_name = "marss_cmd.stan"}
# if (model_type == "Gomp") {cmd_file_name = "marss_gomp_cmd.stan"}
# if (model_type == "MA1") {cmd_file_name = "marss_MA1_cmd.stan"}
# if (model_type == "GP") {cmd_file_name = "marss_gp_cmd_old.stan"}
# if (model_type == "GP_MA1") {cmd_file_name = "marss_gp_ma1_cmd_old.stan"}
# if (model_type == "Q_unconstrained") {cmd_file_name = "marss_Qunconstrained_cmd.stan"}


#=== Run model
file <- file.path(cmdstan_path(), "pinnipeds", cmd_file_name) #  THIS LINE SHOULD CHANGE TO YOUR CMDSTAN DIRECTORY
mod <- cmdstan_model(file)
fit2 <- mod$sample(
  data = data_list_tmp$data,
  seed = 124,
  chains = mcmc_list$n_chain,
  parallel_chains = mcmc_list$n_chain,
  iter_warmup = mcmc_list$n_burn,
  iter_sampling = mcmc_list$n_mcmc,
  adapt_delta=0.8,
  step_size=0.9,
  thin = 1,
  refresh = 100 # print update every 500 iters
)

preds = fit2$summary(variables = "pred")

preds_matrix = cmdstanr_preds_sd(preds,species_i_df)
# saveRDS(preds_matrix,"seabird_preds_matrix.RDS")

# save model diags
diag_summ = fit2$diagnostic_summary()
divergent_perc = sum(diag_summ$num_divergent)/(mcmc_list$n_chain*mcmc_list$n_mcmc)
fit2_sum = fit2$summary()

eff_sample_size = fit2_sum$ess_tail/(mcmc_list$n_chain*mcmc_list$n_mcmc)
n_eff0.01_tmp = sum(eff_sample_size<0.01,na.rm=T)/length(eff_sample_size)

rhat1.1_tmp = sum(fit2_sum$rhat>1.1,na.rm=T)/length(fit2_sum$rhat)

loo_tmp = loo(fit2$draws("log_lik"))

# # Total summed abundance - only coded in for a few cmdstan files, i was lazy
tot_sum_x = fit2$draws("tot_sum_x",format = "matrix")
tot_sum = data.frame(tot.mean = apply(tot_sum_x, 2, mean),
                     tot.sd = apply(tot_sum_x, 2, sd),
                     tot.q025 = apply(tot_sum_x,2,quantile, probs = c(0.25)),
                     tot.q50 = apply(tot_sum_x,2,quantile, probs = c(0.50)),
                     tot.q975 = apply(tot_sum_x,2,quantile, probs = c(0.75)),
                     var="total_abundance_norm") %>%
  mutate(year = as.numeric(colnames(species_i_df)))
# saveRDS(tot_sum,"tot_sum.RDS")

# # LOOIC for model selection
# loo_list_name = paste(species_i,model_type,sep="-")
# loo_all[[loo_list_name]] = loo_tmp
# saveRDS(loo_all,loo_file_name)

model_diags = data.frame(species = "Frigates", #rmu = rmu_i,proc_state = Q, 
                         rhat1.1=rhat1.1_tmp,n_eff0.01=n_eff0.01_tmp,div_p=divergent_perc,
                         loo = loo_tmp$estimates["looic","Estimate"],
                         loo_se = loo_tmp$estimates["looic","SE"]) %>% mutate(model_type = model_type)

# Save Process/Obs Sigmas
if (model_type == "Ueq") {proc_variable_names = c("sigma_process_real")}
if (model_type == "Gomp") {proc_variable_names = c("sigma_process_real","K")}
if (model_type == "MA1") {proc_variable_names = c("sigma_process_real","phi_real")}
if (model_type == "GP") {proc_variable_names = c("sigma_process_real","alpha_real","rho_gp")}
if (model_type == "GP_MA1") {proc_variable_names = c("phi_real","sigma_process_real","alpha_real","rho_gp")}
if (model_type == "Q_unconstrained") {proc_variable_names = c("sigma_process_Q")}


sigma_proc_draws <- fit2$draws(proc_variable_names,format = "matrix")
sigma_proc_tmp = sigma_proc_draws %>% data.frame() %>% pivot_longer(everything(), names_to = "units") %>%
  group_by(units) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            q025 = quantile(value, probs = c(0.025)),
            q50 = quantile(value, probs = c(0.5)),
            q975 = quantile(value, probs = c(0.975))) %>%
  mutate(var="sigma_proc")

# whyy = fit2$draws("Cov",format = "matrix")
sigma_obs_draws <- fit2$draws("sigma_obs",format = "matrix")
sigma_obs_tmp = data.frame(units = levels(fct_inorder(method_vec)), #"sigma_obs", #unique_units_vec, #levels(fct_inorder(method_vec))
                           mean = apply(sigma_obs_draws, 2, mean),
                           sd = apply(sigma_obs_draws, 2, sd),
                           q025 = apply(sigma_obs_draws,2,quantile, probs = c(0.025)),
                           q50 = quantile(sigma_obs_draws, probs = c(0.5)),
                           q975 = apply(sigma_obs_draws,2,quantile, probs = c(0.975)),
                           var="sigma_obs")


error_params_tmp = rbind(sigma_proc_tmp,sigma_obs_tmp) %>% mutate(model_type = model_type)
error_params_tmp$species = species_i 
# error_params = rbind(error_params,error_params_tmp)


# join raw data to plot 
raw_data = species_i_df %>% data.frame() %>% 
  set_colnames(colnames(species_i_df)) %>%
  rownames_to_column("ts") %>% 
  pivot_longer(!ts,
               names_to="year", 
               values_to = "data") %>% 
  mutate(year=as.numeric(year)) 


# extract PRED  from commandstanr
preds = fit2$summary(variables = "pred")
preds_matrix = cmdstanr_preds_sd(preds,species_i_df)
pred_data_tmp = preds_matrix %>% left_join(raw_data) %>% mutate(species = species_i, model_type = model_type) # %>%
  # left_join(tot_sum)

# Plot 
#=== fix the errors to be constant, find the distance (# of years) to the nearest data value
best_pass2_model_fits_dist = pred_data_tmp %>% 
  group_by(ts) %>%
  mutate(distance = NA)

ts_fits_distance = NULL
ts_vec_loop = unique(best_pass2_model_fits_dist$ts)
for (i in 1:length(ts_vec_loop)){
  ts_tmp = ts_vec_loop[i]
  ts_fits_tmp = best_pass2_model_fits_dist %>% filter(ts == ts_tmp)
  
  a <- which(!is.na(ts_fits_tmp$data)) 
  b <- which(is.na(ts_fits_tmp$data))
  
  for (i in 1:length(b)){
    dist_tmp = min(abs(a - b[i]),na.rm=T)
    ts_fits_tmp$distance[b[i]] = dist_tmp
  }
  
  ts_fits_distance = rbind(ts_fits_distance,ts_fits_tmp)
}

# fix the ballooning errors
best_pass2_model_fits_error_fix = ts_fits_distance %>% ungroup() %>% # 
  mutate(sd_tmp = sd) %>%
  mutate(sd_tmp = case_when(distance > 5 ~ NA,
                            .default = sd_tmp)) %>%
  dplyr::group_by(ts) %>%
  # this linearly interpolates NA values INBETWEEN data years, na.rm keeps the ends NA
  # mutate(lo_diff_approx = zoo::na.approx(lo_diff, na.rm=FALSE)) %>% 
  # mutate(hi_diff_approx = zoo::na.approx(hi_diff, na.rm=FALSE)) %>%
  # toggle below on if we want to extend the TS to all years all time
  # ungroup() %>%
  # complete(ts,year, # fill in missing years with NA
  #          fill = list(hi_diff_approx = NA)) %>%
  # group_by(ts) %>%
  # this fills in the ends
  # fill(lo_diff_approx) %>% 
  fill(c(sd_tmp), .direction = "downup") #change to up only if necessary

# plot each time series
p = best_pass2_model_fits_error_fix %>% ggplot() + 
  geom_ribbon(aes(x=year,ymin = mean-1.96*sd_tmp, ymax = mean+1.96*sd_tmp), alpha=0.2) +
  # geom_line(aes(x=year,y=state_fit)) +
  geom_line(aes(x=year,y=mean)) + 
  geom_point(aes(x=year,y=data),color="red")+
  facet_wrap(~ts,scales="free") + theme_bw()
print(p)

# plot total abundance in REAL space
p = tot_sum %>% ggplot() + 
  geom_ribbon(aes(x=year,ymin = tot.q025, ymax = tot.q975), alpha=0.2) +
  geom_line(aes(x=year,y=tot.q50)) +
  theme_bw() + 
  xlim(c(1950,2010)) + ylim(c(0,8e+8))
print(p)
