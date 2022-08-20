library(mcmc) # metrop()

load(url("http://www.stat.umn.edu/geyer/8054/data/nfl2019.rda"))

dim(nfl2019) # 267 games in total
with(nfl2019, sum(week <= 17)) # 256 regular-season game and 11 post-season game
with(nfl2019, sum(road.score == home.score)) # 1 game that has a tie
with(nfl2019, sum(neutral == 1)) # 6 games on neutral sites
with(nfl2019, length(unique(union(home.team, road.team)))) # 32 teams in total

# names for beta, 34 = 32 + 2 in total
varname <- c("ties", "homeadvatange", levels(as.factor(nfl2019$home.team)))
K <- length(varname)
N <- dim(nfl2019)[1]

# Preprocess data
mydata <- nfl2019
# enumerate teams. can convert directly since levels and their orders are same
mydata$home.team <- 2 + as.numeric(as.factor(mydata$home.team))
mydata$road.team <- 2 + as.numeric(as.factor(mydata$road.team))
# regular session or not
mydata$regular <- with(mydata, as.numeric(week <= 17))
# response vector. (homewin, roadwin, tie)
mydata$homewin <- with(mydata, as.numeric(home.score > road.score))
mydata$roadwin <- with(mydata, as.numeric(home.score < road.score))
mydata$tie <- with(mydata, as.numeric(home.score == road.score))
# remove extra information: week, seed, conf
rmindex <- which(colnames(mydata) %in% c("week", "road.seed", "home.seed",
                         "conf", "home.score", "road.score"))
mydata <- mydata[, -rmindex]
# switch postions of two team columns
mydata <- cbind(home.team = mydata$home.team, road.team = mydata$road.team,
                        mydata[, 3:7])
# check with: head(mydata)


# Design model tensor M \in \R^{n \times 3 \times K}
# where n is the number of games, and K is (2 + number of teams)
# (input) gameinfo: a dataframe (home.team, road.team, neutral, regular).
# (cont.) additional columns won't matter as input.
# (cont.) home.team, road.team: number index of team, both should be >= 3
# (cont.) neutral = 1: neutral site; regular = 1; regular session game
# (output) M: modeltensor as described above
modeltensor <- function(gameinfo) {
    N <- dim(gameinfo)[1]
    M <- array(0, c(3, K, N))

    for (n in 1:N) {
        home.team <- gameinfo[n, 1]
        road.team <- gameinfo[n, 2]
        neutral <- gameinfo[n, 3]
        regular <- gameinfo[n, 4]

        result <- M[, , n]
        # effect of tie
        result[3, 1] <- 1
        # effect of home advantage
        if (neutral == 0)
            result[c(1, 3), 2] <- c(1, 1 / 2)
        # effect teams
        result[c(1, 3), home.team] <- c(1, 1 / 2)
        result[c(2, 3), road.team] <- c(1, 1 / 2)
        # correction of post-session game (I think use -Inf should be better)
        if (regular == 0)
            result[3, ] <- NA # use 0 for canonical stat and -Inf for \theta

        M[, , n] <- result
    }
    return(M)
}

# model tensor for our dataset: M \in \R^{n \times 3 \times K}
# (cont.) with M[, , i] the model matrix for the ith game
M <- modeltensor(mydata) # original M with NA
M_0 <- M
M_0[which(is.na(M_0))] <- 0 # for calculating cano stat


# unnormalized log density with respect to \beta, that is,
# \ell (\beta) = \innerproduct{M^\T y}{\beta} - \sum_{i=1}^n c_{sub}(\beta)
# (input) beta: a numerical vector of size K, with \beta as latent parameters
log.unnormalized.density <- function(beta) {
    stopifnot(is.numeric(beta))
    stopifnot(is.finite(beta))

    # canonical statistics
    canostat <- NULL
    for (k in 1:K) {
        temp <- sum(diag(M_0[, k, ] %*% as.matrix(mydata[, 5:7])))
        canostat <- c(canostat, temp)
    }

    # cumulant part
    cumulant <- 0
    for (n in 1:N) {
        theta_i <- M[, , n] %*% beta

        # deal with non-regular session
        if (mydata[n, "regular"] == 0)
           theta_i[3, 1] <- -Inf # always -Inf so 0 prob for tie in post-session

        cumulant_i <- log(sum(exp(theta_i)))
        cumulant <- cumulant + cumulant_i
    }

    # combine them together
    sum(canostat * beta) - cumulant
}

# unnormalized log prior
log.unnormalized.prior <- function(beta) {
    stopifnot(is.numeric(beta))
    stopifnot(is.finite(beta))

    logp <- ifelse(beta < 0, beta - log1p(exp(beta)), - log1p(exp(- beta)))
    logq <- ifelse(beta < 0, - log1p(exp(beta)), - beta - log1p(exp(- beta)))
    sum(logp) + sum(logq)
}
# unnormalized log posterior
log.unnormalized.posterior <- function(beta)
    log.unnormalized.density(beta) + log.unnormalized.prior(beta)


# MCMC
set.seed(1234) # for reproducibility
mout <- metrop(log.unnormalized.posterior, rep(0, K),
    nbatch = 150, blen = 150, scale = 0.21)
mout$accept

mout <- metrop(mout, scale = 0.215)
mout$accept

mout <- metrop(mout, scale = 0.22)
mout$accept
t.test(mout$accept.batch)$conf.int

# save(mout, file = "mout.RData")
# load(file = "mout.Rdata")

plot(ts(mout$batch[, 1]), main = "Batch Means for Tie coefficient")
plot(ts(mout$batch[, 2]), main = "Batch Means for Home Advantage Coefficient")
acf(ts(mout$batch[, 2]),
        main = "Autocorrelation of Home Advantage Coefficient")
plot(ts(mout$batch[, 3]))



# Specific question 1
n <- 267 # game index
outfun <- function(beta) {
    theta_i <- M[, , n] %*% beta
    # deal with non-regular session
    if (mydata[n, "regular"] == 0)
        theta_i[3, 1] <- -Inf # always -Inf so 0 prob for tie in post-session

    exp(theta_i) / sum(exp(theta_i))
}

mout_1 <- metrop(mout, blen = 1000, outfun = outfun)
# save(mout_1, file = "mout_Q1.RData")
# load(file = "mout_Q1.RData")

# prob of homewin, roadwin and tie
colMeans(mout_1$batch)


# Specific question 2
# related teams and their index
finalteams <- with(nfl2019[257:267, ], unique(union(home.team, road.team))) # 12
teamindex <- which(varname %in% finalteams)
teamname <- varname[teamindex]

teamseed <- rep(0, 12)
teamconf <- rep("AFC", 12)
temp <- nfl2019[257:266, ] # no need for the last game with two diff confs
for (i in 1:dim(temp)[1]) {
    index <- which(teamname == temp$road.team[i])
    teamseed[index] <- temp$road.seed[i]
    teamconf[index] <- temp$conf[i]

    index <- which(teamname == temp$home.team[i])
    teamseed[index] <- temp$home.seed[i]
    teamconf[index] <- temp$conf[i]
}
# aggregate information of teams in the playoffs
finalgames <- data.frame(teamname, teamindex, teamseed, teamconf)


# simulate the playoffs recursively
outfun2 <- function(beta) {
    # within AFC
    ## result of week 18
    w18_AFC1 <- ifelse(beta[2] + beta[32] > beta[6], 32, 6)
    w18_AFC2 <- ifelse(beta[2] + beta[24] > beta[33], 24, 33)
    ### record seed information
    seed1 <- with(finalgames, teamseed[teamindex == w18_AFC1])
    seed2 <- with(finalgames, teamseed[teamindex == w18_AFC2])
    ## result of week 19
    if (seed1 < seed2) {
        w19_AFC1 <- ifelse(beta[2] + beta[27] > beta[w18_AFC1], 27, w18_AFC1)
        w19_AFC2 <- ifelse(beta[2] + beta[12] > beta[w18_AFC2], 12, w18_AFC2)
    }
    if (seed1 > seed2) {
        w19_AFC1 <- ifelse(beta[2] + beta[27] > beta[w18_AFC2], 27, w18_AFC2)
        w19_AFC2 <- ifelse(beta[2] + beta[12] > beta[w18_AFC1], 12, w18_AFC1)
    }
    ### record seed information
    seed1 <- with(finalgames, teamseed[teamindex == w19_AFC1])
    seed2 <- with(finalgames, teamseed[teamindex == w19_AFC2])
    ## result of week 20
    if (seed1 < seed2)
        w20_AFC <- ifelse(beta[2] + beta[w19_AFC1] > beta[w19_AFC2],
                          w19_AFC1, w19_AFC2)
    if (seed1 > seed2)
        w20_AFC <- ifelse(beta[2] + beta[w19_AFC2] > beta[w19_AFC1],
                          w19_AFC2, w19_AFC1)


    # within NFC
    ## result of week 18
    w18_NFC1 <- ifelse(beta[2] + beta[29] > beta[34], 29, 34)
    w18_NFC2 <- ifelse(beta[2] + beta[16] > beta[30], 16, 30)
    ### record seed information
    seed1 <- with(finalgames, teamseed[teamindex == w18_NFC1])
    seed2 <- with(finalgames, teamseed[teamindex == w18_NFC2])
    ## result of week 19
    if (seed1 < seed2) {
        w19_NFC1 <- ifelse(beta[2] + beta[3] > beta[w18_NFC1], 3, w18_NFC1)
        w19_NFC2 <- ifelse(beta[2] + beta[22] > beta[w18_NFC2], 22, w18_NFC2)
    }
    if (seed1 > seed2) {
        w19_NFC1 <- ifelse(beta[2] + beta[3] > beta[w18_NFC2], 3, w18_NFC2)
        w19_NFC2 <- ifelse(beta[2] + beta[22] > beta[w18_NFC1], 22, w18_NFC1)
    }
    ### record seed information
    seed1 <- with(finalgames, teamseed[teamindex == w19_NFC1])
    seed2 <- with(finalgames, teamseed[teamindex == w19_NFC2])
    ## result of week 20
    if (seed1 < seed2)
        w20_NFC <- ifelse(beta[2] + beta[w19_NFC1] > beta[w19_NFC2],
                          w19_NFC1, w19_NFC2)
    if (seed1 > seed2)
        w20_NFC <- ifelse(beta[2] + beta[w19_NFC2] > beta[w19_NFC1],
                          w19_NFC2, w19_NFC1)

    # final game between AFC and NFC! week 21!
    w21 <- ifelse(beta[w20_AFC] > beta[w20_NFC], w20_AFC, w20_NFC)

    result <- rep(0, dim(finalgames)[1])
    index <- with(finalgames, which(teamindex == w21))
    result[index] <- 1
    names(result) <- finalgames$teamname
    result
}

mout_2 <- metrop(mout, blen = 1000, outfun = outfun2)
# save(mout_2, file = "mout_Q2.RData")
# load(file = "mout_Q2.RData")
winprob <- colMeans(mout_2$batch)
names(winprob) <- finalgames$teamname

# take results that are significant
foo <- as.ts(mout_2$batch)
colnames(foo) <- finalgames$teamname
foo <- foo[, colMeans(foo) > 0.1]

# MCSE
apply(foo, 2, sd) / sqrt(mout$nbatch)

# basic diagnostic plots
plot(foo)
acf(foo)

# library(knitr)
# bar <- rbind(colMeans(foo), apply(foo, 2, sd) / sqrt(mout$nbatch))
# rownames(bar) <- c("probability of being best", "MCSE")
# kable(bar, digits = 2, caption = "Estimated Posterior Probability of Being the Best Team in the Conference.  Probabilities to not add to one; there was small probability of other teams being the best.")