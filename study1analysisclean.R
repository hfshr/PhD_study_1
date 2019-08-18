load(file = "datavars/dataprepro.rds")
load(file = "datavars/finalrst.rds")
rm(missingcount, missingid, missingtot)
library(tidyverse)
library(zoo)
library(bnlearn)
library(Rgraphviz)

## function to do median split
splitr <- function(x){cut(x, breaks = c(min(x), median(x), max(x)), labels = c("Low", "High"), 
                          include.lowest = TRUE, ordered_result = TRUE)}
## function to prepare data for network
networkprepare <- function(networkdata){
  filter(networkdata, time == 1) %>% select(-time) -> t1
  filter(networkdata, time == 2) %>% select(-time) -> t2
  filter(networkdata, time == 3) %>% select(-time) -> t3
  explans <- networkdata %>% 
    select(c(id, pi:hours)) %>% 
    names()
  merge(t1, t2, by = explans, suff = c("_1", "_2")) %>% 
    select(-id) -> t12
  merge(t2, t3, by = explans, suff = c("_1", "_2")) %>% 
    select(-id) -> t23
  tallc = rbind(t12, t23)
} # generates tallc files
networkmaker <- function(names){
  networkdata <- datanet %>% 
    select(names) %>% 
    mutate_at(vars(BIS:FFFS, nlebase), splitr)
  tallc <- networkprepare(networkdata)
  
  g1 = tallc %>% select(ends_with("_1")) %>% names()
  g2 = tallc %>% select(ends_with("_2")) %>% names()
  
  explans <- networkdata %>% 
    select(pi:hours, nlebase) %>% 
    names()
  
  bl <-  expand.grid(from = g2, to = g1, stringsAsFactors = FALSE)
  bl <- rbind(bl, expand.grid(from = c(g1, g2), 
                              to = explans, stringsAsFactors = FALSE)) 
  bl <- rbind(bl, data.frame(from = c("clevel", "nlebase"),
                             to = c("gender", "ind_team"), stringsAsFactors = FALSE)) 
  bl <- rbind(bl, data.frame(from = c(g1[2:length(g1)]),
                             to = c(g2[2:length(g2)]), stringsAsFactors = FALSE))
  
  
  wl = data.frame(from = c(paste(vars_select(names(tallc), starts_with("n"))[2]), 
                           paste(vars_select(names(tallc), starts_with("n"))[3])),
                  to = c("injured_1", "injured_2"))
  
  str = boot.strength(tallc, algorithm = "tabu",
                      algorithm.args = list(blacklist = bl, 
                                            #whitelist = wl,
                                            tabu = 50))
  
  
  return(list(g1 = g1,
              g2 = g2,
              explans = explans,
              tallc = tallc,
              str = str))
  
}

## data for network
datanet <- datanew %>% 
  select(-bis, -bas, -fffs, -delta, -post, -pre, -change) %>% 
  left_join(., finalrst, by = c("id", "time")) %>%  # updated 51 item rst
  filter(variable == 1) %>% 
  rename("balance" = gt,
         "pi" = pic ) %>% 
  mutate(injured = factor(ifelse(injured == 0, "healthy", "injured")),
         clevel = factor(clevel, levels = c("club_university_county",
                                            "national", "international")),
         clevel = fct_collapse(clevel, "national_international" = c("national", "international")),
         pi = factor(pi),
         gender = factor(gender),
         sportg = factor(sportg),
         hours = factor(ifelse(hours < median(hours), "Low", "High"), levels = c("Low", "High")),
         sportg = factor(sportg),
         ind_team = fct_collapse(sportg, "individual" = c("athletics", "gym", "other"), 
                                 "team" = c("football", "hockey", "netball", "cricket", 
                                            "rugby", "basketball"))
  ) %>%  
  group_by(id) %>% 
  mutate(nlec = cumsum(nle),
         tlec = cumsum(tle),
         nlebase = ifelse(time == 1, nle, NA),
         nlebase = na.locf(nlebase),
         injnum = ifelse(injured == "injured", 1, 0),
         injcount = sum(injnum),
         oneinj = ifelse(injcount >= 1, 1, 0),
         oneinj = factor(oneinj, levels = c(0 ,1), labels = c("healty", "injured"))) %>% 
  ungroup()  %>% 
  group_by(time) %>% 
  mutate(nlel = log(nlec +1),
         tlel = log(tlec)) %>% 
  ungroup() %>% 
  mutate(negsev = ifelse(negsev == "NaN", 0, negsev)) %>% 
  ungroup() %>% 
  mutate(rmssd = round(log(rmssd), digits = 2)) %>% 
  group_by(time) %>% 
  mutate(totneg_d = splitr(totneg),
         "nlelg" = splitr(log(nlec +1)),
         "tlelg" = splitr(log(tlec)),
         nlel = splitr(log(nle + 1)),
         tlel = splitr(log(tle + 1))) %>% 
  ungroup()
rm(finalrst)

save(datanet, file = "datavars/networkdata.rds")
load("datavars/networkdata.rds")


## insert new explanatory variables between `pic` and `hours` to ensure later functions work
networkdata <- datanet %>% 
  select(id, time, injured, nlelg, pi, nlebase, gender, ind_team, clevel, hours, stiffness,
         RI, BIS, rmssd, balance, FFFS) %>% 
  mutate_at(vars(stiffness:FFFS, nlebase), splitr)

## prepare 2 turn network
tallc <- networkprepare(networkdata)

## table for included variables 
nlesplitr <- function(timenum){
  datanet %>% 
    select(time, nlec) %>%
    filter(time == timenum) %>%
    group_by(time) %>% 
    mutate(lower = log(nlec +1)) %>% 
    mutate_at(vars(lower), 
              function(x) cut(x, breaks = c(min(x), median(x), max(x)), 
                              include.lowest = TRUE, ordered_result = TRUE)) %>% 
    map(~as.numeric(sub('.(.+),.+', '\\1', levels(.x)))[2]) %>% 
    as_tibble() %>% 
    select(-nlec)  %>% 
    mutate(var = paste("nlelg_", timenum, sep = "")) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(nlec) %>% 
            mutate(nlec = log(nlec +1)) %>% 
            map(~min(.)) %>% 
            as_tibble() %>% 
            gather(varr, min) %>% 
            select(-varr)) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(nlec) %>% 
            mutate(nlec = log(nlec +1)) %>% 
            map(~max(.)) %>% 
            as_tibble() %>% 
            gather(varr, max) %>% 
            select(-varr)) %>% 
    mutate(max = round(max, digits = 2)) %>% 
    mutate(Low = paste(min, lower, sep = "-"),
           High = paste(">", lower, sep = ""),
           High = paste(High, max, sep = "-")) %>% 
    mutate(Definition = "") %>% 
    select(var, Definition, Low, High) %>% 
    rename(state_1 = "Low", 
           state_2 = "High")
  
}
tlesplitr <-  function(timenum){
  datanet %>% 
    select(time, tlec) %>%
    filter(time == timenum) %>%
    group_by(time) %>% 
    mutate(lower = log(tlec)) %>% 
    mutate_at(vars(lower), 
              function(x) cut(x, breaks = c(min(x), median(x), max(x)), 
                              include.lowest = TRUE, ordered_result = TRUE)) %>% 
    map(~as.numeric(sub('.(.+),.+', '\\1', levels(.x)))[2]) %>% 
    as_tibble() %>% 
    select(-tlec)  %>% 
    mutate(var = paste("tlelg_", timenum, sep = "")) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(tlec) %>% 
            mutate(tlec = log(tlec)) %>% 
            map(~min(.)) %>% 
            as_tibble() %>% 
            gather(varr, min) %>% 
            select(-varr)) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(tlec) %>% 
            mutate(tlec = log(tlec)) %>% 
            map(~max(.)) %>% 
            as_tibble() %>% 
            gather(varr, max) %>% 
            select(-varr)) %>% 
    mutate(max = round(max, digits = 2),
           min = round(min, digits = 2)) %>% 
    mutate(Low = paste(min, lower, sep = "-"),
           High = paste(">", lower, sep = ""),
           High = paste(High, max, sep = "-")) %>% 
    mutate(Definition = "") %>% 
    select(var, Definition, Low, High) %>% 
    rename(state_1 = "Low", 
           state_2 = "High") 
} 

names <- datanet %>% 
  select(nlebase, stiffness:FFFS) %>% 
  names()

varcutoffs <- tallc %>% 
  select(pi:hours)%>%
  select(-nlebase) %>% 
  map(~levels(.)) %>% 
  as_tibble() %>% 
  gather(var, state) %>% 
  mutate(ind = rep(c("state_1", "state_2"), length.out = n())) %>% 
  spread(ind, state)%>% 
  mutate(Definition = c("Current competitive level",
                        "Gender of the participant",
                        "Number of hours spent training per week",
                        "Participate in an individual or team based sport",
                        "Previous injury - Whether an injury had been sustained in the previous 12 months prior to the study")) %>% 
  mutate_at(vars(state_1, state_2), funs(str_to_title(.))) %>% 
  select(var, Definition, state_1, state_2) %>% 
  rbind(., datanet %>% 
          select(id, time, injured, pi, nlebase, gender, ind_team, clevel, hours, stiffness, 
                 RI, RR, I, GDP, BIS, rmssd, sdnn, bal_asym, balance, FFFS) %>% 
          mutate_at(vars(stiffness:FFFS, nlebase), 
                    function(x) cut(x, breaks = c(min(x), median(x), max(x)), 
                                    include.lowest = TRUE, ordered_result = TRUE)) %>% 
          select(stiffness:FFFS, nlebase) %>% 
          map(~as.numeric(sub('.(.+),.+', '\\1', levels(.x)))[2]) %>% 
          as_tibble()%>%  
          gather(var, lower) %>% 
          {. ->> namess} %>% 
          cbind(., datanet %>% 
                  select(namess$var) %>% 
                  map(~min(.)) %>% 
                  as_tibble() %>% 
                  gather(varr, min) %>% 
                  select(-varr)) %>% 
          cbind(., datanet %>% 
                  select(namess$var) %>% 
                  map(~max(.)) %>% 
                  as_tibble() %>% 
                  gather(varr, max) %>% 
                  select(-varr)) %>% 
          mutate(Low = paste(min, lower, sep = "-"),
                 High = paste(">", lower, sep = ""),
                 High = paste(High, max, sep = "-")) %>% 
          select(var, Low, High) %>% 
          mutate(var = factor(var, levels = c("nlebase", "FFFS", "BIS", "RI",
                                              "RR", "I", "GDP", "stiffness", "rmssd",
                                              "sdnn", "bal_asym", "balance"))) %>% 
          arrange(var) %>% 
          mutate(Definition = c("Untransformed NLE at TP 1",
                                "Fight-Flight-Freeze System",
                                "Behavioural Inhibition System",
                                "Reward Interest",
                                "Reward reactivity",
                                "Impulsivity",
                                "Goal drive persistence",
                                "Sum of all stiffness locations",
                                "Root mean squared difference of successive RR intervals",
                                "Standard deviation of RR series",
                                "Percentage difference between left and right leg balance score",
                                "Total balance score")) %>% 
          rename(state_1 = "Low", 
                 state_2 = "High") %>% 
          select(var, Definition, state_1, state_2))%>% 
  rbind(., map_df(1:3, nlesplitr)) %>% 
  rbind(., map_df(1:3, tlesplitr)) %>% 
  rowid_to_column("tempid") %>% 
  mutate(state_1 = ifelse(tempid > 5, paste(state_1, "(Low)"), state_1),
         state_2 = ifelse(tempid > 5, paste(state_2, "(High)"), state_2),
         state_1 = ifelse(tempid == 3 ,"0-9 (Low)", state_1),
         state_2 = ifelse(tempid == 3, ">9-35 (High)", state_2)) %>% 
  select(-tempid)

varcutoffs[18,2] <- "Log NLE at TP 1"
varcutoffs[19,2] <- "Log NLE at TP 2"
varcutoffs[20,2] <- "Log NLE at TP 3"
varcutoffs[21,2] <- "Log TLE at TP 1"
varcutoffs[22,2] <- "Log TLE at TP 2"
varcutoffs[23,2] <- "Log TLE at TP 3"



## The network - extract names for blacklist + plotting

g1 <-  tallc %>% select(ends_with("_1")) %>% names()
g2 <-  tallc %>% select(ends_with("_2")) %>% names()
  
  explans <- networkdata %>% 
    select(pi:hours, nlebase) %>% 
    names()
  
  bl <-  expand.grid(from = g2, to = g1, stringsAsFactors = FALSE) %>% 
    rbind(., expand.grid(from = c(g1, g2), to = explans, stringsAsFactors = FALSE)) %>% 
    rbind(., data.frame(from = c("clevel", "nlebase", "nlebase"),
                            to = c("gender", "ind_team", "gender"), stringsAsFactors = FALSE)) %>% 
    rbind(., data.frame(from = c(g1[2:length(g1)]),
                            to = c(g2[2:length(g2)]), stringsAsFactors = FALSE))

wl = data.frame(from = c("nlelg_1", "nlelg_2"),
                to = c("injured_1", "injured_2"))


library(parallel)
cl = makeCluster(2, type = "SOCK")

set.seed(125)
str = boot.strength(tallc, algorithm = "tabu",
                    algorithm.args = list(blacklist = bl,
                                          whitelist = wl,
                                          tabu = 500),
                    cluster = cl)

avg30 = averaged.network(str, threshold = 0.30)
fitted = bn.fit(avg30, tallc)

## plot the big network
gR = strength.plot(avg30, str, 
                   shape = "rectangle", 
                   render = FALSE, 
                   groups = list(g1, g2, explans),
                   layout = "dot")

renderGraph(gR)


nodeRenderInfo(gR)$fill[g1] = "deepskyblue"
nodeRenderInfo(gR)$fill[g2] = "tomato"
nodeRenderInfo(gR)$fill["injured_1"] = "gold"
nodeRenderInfo(gR)$fill["injured_2"] = "gold"
gR <- renderGraph(gR)

save(varcutoffs, avg30, datanet, fitted, gR, networkdata, str, tallc, file = "datavars/tempnetworks.rds")
load("datavars/tempnetworks.rds")


### temp test for non cumlative nle

mb(avg30, "injured_1")
var <- c(mb(avg30, "injured_1"))
query1 <- map_dfr(var, ~newprobtable(tallc,
                                     .x,
                                     "injured_1",
                                     "injured",
                                     avg30)) %>% 
  gather(variable, state, -prob, - outcome) %>% 
  drop_na() %>% 
  select(variable, state, prob) %>% 
  spread(state, prob) %>% 
  select(variable, Low, High) %>% 
  dendroTools::round_df(2)






## compare netwrok scores for different models 

nlename <- list(nlelg = datanet %>% 
                  select(id, time, injured, nlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                         RR, I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names(),
                tlelg = datanet %>% 
                  select(id, time, injured, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                         RR, I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names(),
                both = datanet %>% 
                  select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                         RR, I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names()
)
hrvnames <- list(rmssd = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RI, I, GDP, RI, rmssd, balance, FFFS) %>% names(),
                 sdnn = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RI,I, GDP, RI, sdnn, balance, FFFS) %>% names(),
                 both = datanet %>% 
                   select(id, time, injured, nlelg,  tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RI,I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names())

rstnames <- list(BAS = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          BAS, rmssd,sdnn, balance, FFFS) %>% names(),
                 RI = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RI, rmssd,sdnn, balance, FFFS) %>% names(),
                 GDP = datanet %>% 
                   select(id, time, injured,nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          GDP, rmssd,sdnn, balance, FFFS) %>% names(),
                 I = datanet %>% 
                   select(id, time, injured, nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          I, rmssd,sdnn, balance, FFFS) %>% names(),
                 RR = datanet %>% 
                   select(id, time, injured,nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RR, rmssd,sdnn, balance, FFFS) %>% names(),
                 ALL = datanet %>% 
                   select(id, time, injured, nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RR, I, GDP, RI, rmssd,sdnn, balance, FFFS) %>% names())

balancenames <- list(both = datanet %>% 
                       select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                              RI, I, GDP, RI, rmssd, sdnn, balance, bal_asym, FFFS) %>% names(),
                     balance = datanet %>% 
                       select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                              RI,I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names(),
                     asym = datanet %>% 
                       select(id, time, injured, nlelg,  tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                              RI,I, GDP, RI, sdnn, rmssd, bal_asym, FFFS) %>% names())

stiffness <- list(all = datanet %>% 
                    select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                           frequency, decrement, RI, I, GDP, RI, rmssd, sdnn, balance, FFFS) %>% names(),
                  stiffness = datanet %>% 
                    select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                           RI,I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names(),
                  frequency = datanet %>% 
                    select(id, time, injured, nlelg,  tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, frequency, 
                           RI,I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names(),
                  elasticity = datanet %>% 
                    select(id, time, injured, nlelg,  tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, decrement, 
                           RI,I, GDP, RI, sdnn, rmssd, balance, FFFS) %>% names())



library(tidyselect)

scorer <- function(names){
  test <- map(names, ~networkmaker(.x))
  avg <- map(.x = names(test), ~averaged.network(test[[.x]][[5]], threshold = 0.35)) %>% 
    set_names(c(names(test)))
  "Variable" <- names(test)
  x <- map_df(names(test), ~data.frame(score = score(avg[[.x]], test[[.x]][[4]], type = "bic"))) %>% 
    cbind(Variable, . ) %>% 
    arrange(-score)
}

le <- scorer(nlename)
rstmodel <- scorer(rstnames)
hrv <- scorer(hrvnames)
balance <- scorer(balancenames)

#myoton <- scorer(stiffness)


save(le, rstmodel, hrv, balance, file = "datavars/modeleval.rds")

###############################################################################################
                              #Old code but maybe still useful#
###############################################################################################

## network queries
# function to use for combinations of states with multiple variables
newprobtable <- function(tallc, vars, outcome, state, model) {
  all.levels = if(any(length(vars) > 1)) { 
    lapply(tallc[, (names(tallc) %in% vars)], levels) 
  } else {
    all.levels <- tallc %>% 
      select(vars) %>% 
      sapply(levels) %>% 
      as_tibble() 
  } # makes the code work for when only one variable is used as evidence
  combos <- do.call("expand.grid", c(all.levels, list(stringsAsFactors = FALSE)))  # al combiations
  
  str1 <- "" 
  for (i in seq(nrow(combos))) {
    str1[i] = paste(combos %>% names(), " = '",
                    combos[i, ] %>% sapply(as.character), "'",
                    sep = "", collapse = ", ")
  } # generate character strings for all combinations
  str1 <- rep(str1, times = length(outcome)) # repeat the string for more than one outcome
  str1 <- paste("list(", str1, ")", sep = "")
  
  all.levels.outcome = if(any(length(outcome) > 1)) {
    lapply(tallc[, (names(tallc) %in% outcome)], levels)
  } else {
    all.levels <- tallc %>% 
      select(outcome) %>% 
      sapply(levels) %>% 
      as_tibble()
  } # repeat loop for outcome variables (can have more than one outcome)
  
  combos.outcome <- do.call("expand.grid", c(all.levels.outcome, list(stringsAsFactors = FALSE)))
  
  str3 = rep(paste("(", outcome, " == '", state, "')", sep = ""), each = length(str1)/length(outcome))  # repeat each outcome for the length of combos
  
  fitted <-  bn.fit(avg30, tallc, method = "bayes", iss = 1) # fit the model with bayes method
  cmd = paste("cpquery(fitted, ", str3, ", ", str1, ", method = 'lw', n = 4000)", sep = "") # join all elements of string together
  prob <-  rep(0, length(str1)) # empty vector for probabilities 
  for (i in seq(length(cmd))){
    prob[i] <- eval(parse(text = cmd[i]))
  } # for each combination of strings, what is the probability of outcome
  test <- cbind(combos, prob) %>% 
    mutate(outcome = str3)
  # output
  # predict <- lm(prob ~ ., data = combos)
  # return_list <- list("combos" = combos, "result" = test, "model" = predict,
  #                     "cmd" = cmd)
  # return(return_list)
}

var <- c("balance_1", "stiffness_1", "nlelg_1", "hours")
query1 <- map_dfr(var, ~newprobtable(tallc,
                                     .x,
                                     "injured_1",
                                     "injured",
                                     avg30)) %>% 
  gather(variable, state, -prob, - outcome) %>% 
  drop_na() %>% 
  select(variable, state, prob) %>% 
  spread(state, prob) %>% 
  select(variable, Low, High) %>% 
  dendroTools::round_df(2)


toptail <- function(data){
  data %>% 
    top_n(3, prob) %>% 
    rbind(., data %>% 
            top_n(3, -prob))
} #function to get hightest and lowest values (n = 3)

query2 <- newprobtable(tallc, 
                                 c("balance_1", "stiffness_1", "nlelg_1", "hours"),
                                 "injured_1",
                                 "injured",
                                 avg30
) %>% 
  toptail() %>% 
  arrange(-prob) %>% 
  select(outcome, everything()) %>% 
  dendroTools::round_df(2)

query2a <- query2[c(1,3),] %>% 
  dendroTools::round_df(2)
  
var2 <- c("balance_2", "stiffness_2", "nlelg_2", "FFFS_1")
query3 <- map_dfr(var2, ~newprobtable(tallc,
                                               .x,
                                               "injured_2",
                                               "injured",
                                               avg30)) %>% 
  gather(variable, state, -prob, - outcome) %>% 
  drop_na() %>% 
  select(variable, state, prob) %>% 
  spread(state, prob) %>% 
  select(variable, Low, High) %>% 
  dendroTools::round_df(2)

query4 <- newprobtable(tallc, 
                       var2,
                       "injured_2",
                       "injured",
                       avg30
) %>% 
  toptail() %>% 
  arrange(-prob) %>% 
  select(outcome, everything()) %>% 
  dendroTools::round_df(2)



map(mb(avg30, "injured_1"), 
        ~newprobtable(tallc, .x, "injured_1", "injured", avg30))


comp <- newprobtable(tallc, vars = names(tallc[c(1:6, 8:14)]), outcome = "injured_1",
                     state = "injured") 
injured1comp <- comp %>% 
  top_n(1, prob) %>% 
  as_tibble() %>% 
  mutate(id = "highrisk",
         prob = round(prob, digits = 8)) %>% 
  rbind(., comp %>% 
          top_n(1, -prob) %>% 
          as_tibble() %>% 
          mutate(id = "lowrisk",
                 prob = round(prob, digits = 8))) %>% 
  mutate(prob = round(prob, digits = 2)) %>% 
  gather(variable, value, -id) %>% 
  spread(id, value) %>% 
  mutate(variable = factor(variable, levels = c("outcome", "prob", "pi", "clevel", "gender", "hours", 
                                      "ind_team", "nlebase", 
                                      "stiffness_1", "balance_1", "RI_1", "BIS_1", "FFFS_1", 
                                      "nlelg_1",  "rmssd_1"))) %>% 
  arrange(variable) 

save(comp, injured1comp, file = "datavars/largequeries.rds")


names(tallc)
comp2 <- newprobtable(tallc, vars = names(tallc[c(1:14, 16:22)]), outcome = "injured_2",
                              state = "injured") 



## test the probability that status does not change between _1 and _2
library(tidyselect)

fitted <- bn.fit(avg30, tallc)

g1new <- paste(g1[2:8], "== 'High'", sep = " ")
g2new <- paste(g2[2:8], "== 'High'", sep = " ")
cmd = paste("cpquery(fitted, ", "event = (", g2new,  "), ", "evidence = (", g1new, "))", sep = "") 

x <- map_df(cmd, ~tibble(name =  eval(parse(text = .x)))) %>% 
  cbind(., as_vector(g1[2:8]))


one <- tallc %>% 
  select(vars_select(names(tallc), ends_with("_1"))) %>% 
  select(-1) %>% 
  rownames_to_column("id") %>% 
  gather(var_1, state_1, -id) %>% 
  cbind(., tallc %>% 
          select(vars_select(names(tallc), ends_with("_2"))) %>% 
          select(-1) %>% 
          rownames_to_column("id") %>% 
          gather(var_2, state_2, -id) %>% 
          select(-id)) %>%  
  mutate(bigvar = str_remove(var_1, "_1")) 
  
counter <- function(state1, state2){
  one %>% 
    group_by(bigvar) %>% 
    count(state_1 == state1 & state_2 == state2) %>% 
    `colnames<-`(c("var", "case", "count")) %>% 
    filter(case == TRUE) %>% 
    ungroup() %>% 
    mutate(perc = count/nrow(tallc)*100,
           state = paste(state1, state2, sep = "/"))
}

combo <- counter("High", "Low") %>% 
  rbind(., counter("Low", "High")) %>% 
  select(-case) %>% 
  mutate(perc = round(perc, digits = 2)) %>% 
  arrange(var)

comboequal <- counter("High", "High") %>% 
  rbind(., counter("Low", "Low")) %>% 
  select(-case) %>% 
  mutate(perc = round(perc, digits = 2)) %>% 
  arrange(var)


networkmakerbl <- function(names, bl){
  networkdata <- datanet %>% 
    select(names) %>% 
    mutate_at(vars(stiffness:FFFS, nlebase), splitr)
  
  tallc <- networkprepare(networkdata)
  
  g1 = tallc %>% select(ends_with("_1")) %>% names()
  g2 = tallc %>% select(ends_with("_2")) %>% names()
  
  explans <- networkdata %>% 
    select(pi:hours, nlebase) %>% 
    names()

  str = boot.strength(tallc, algorithm = "tabu",
                      algorithm.args = list(blacklist = bl, whitelist = wl, tabu = 50))
  
  return(list(g1,g2, explans, tallc, str))
  
}


blist <- list(bl_nojoins = expand.grid(from = g2, to = g1, stringsAsFactors = FALSE) %>% 
                rbind(., expand.grid(from = c(g1, g2), to = explans, stringsAsFactors = FALSE)) %>% 
                rbind(., data.frame(from = c("clevel", "nlebase", "nlebase"),
                                    to = c("gender", "ind_team", "gender"), stringsAsFactors = FALSE)) %>% 
                rbind(., data.frame(from = c(g1[2:length(g1)]),
                                    to = c(g2[2:length(g2)]), stringsAsFactors = FALSE)),
              
              bl_alljoins = expand.grid(from = g2, to = g1, stringsAsFactors = FALSE) %>% 
                rbind(., expand.grid(from = c(g1, g2), to = explans, stringsAsFactors = FALSE)) %>% 
                rbind(., data.frame(from = c("clevel", "nlebase", "nlebase", "nlebase", "nlebase"),
                                    to = c("gender", "ind_team", "gender", "nlelg_1", "nlelg_2"), stringsAsFactors = FALSE)))

test <- map(.x = blist, ~networkmakerbl(networkdata %>% names(), .x))
avg <- map(.x = c(1:length(test)), ~averaged.network(test[[.x]][[5]], threshold = 0.35))
map(c(1:length(test)), ~score(avg[[.x]], test[[.x]][[4]], type = "bic")) %>% 
  set_names(names(blist))



## changes


