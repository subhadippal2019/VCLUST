
fun<-function(x){exp(-(4-x^2)^2)}
dMinus_Log_Fun<-function(x){ -(2 * (2 * x * (4 - x^2)))}
d2Minus_Log_Fun<-function(x){-(2 * (2 * (4 - x^2) - 2 * x * (2 * x)))}

x <- rMARS(100,"exp(-(4-x^2)^2)",-Inf,Inf, c(-2.5,0,2.5),c(-2/sqrt(3),2/sqrt(3)))


x1 <- rMARSs(100,fun=fun,dMinus_Log_Fun=dMinus_Log_Fun, d2Minus_Log_Fun=d2Minus_Log_Fun ,-Inf,Inf, c(-2.5,0,2.5),c(-2/sqrt(3),2/sqrt(3)))

rMARSs<-function (n, fun,dMinus_Log_Fun,d2Minus_Log_Fun , min = -Inf, max = Inf, sp, infp, m = 10^(-4))
{
  browser()
  sp <- sort(sp)
  infp <- sort(infp)
  #if (!is.character(formula))
    #stop("Density function is inappropriate, please look up examples for help")
  if (n <= 0)
    stop("Length of sequence shouble be larger than 0")
  if (min >= max)
    stop("Minimum of domain shouble be larger than maximum")
  ltuInf <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(pandl[x, 2][[1]], infp[1])
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(pandl[x, 3][[1]])
    crv <- numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      A = matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2,
                 byrow = T)
      b = -c(int[i], int[i + 1])
      crp[i] = solve(A, b)[1]
      crv[i] = solve(A, b)[2]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  ltuFin <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(min, pandl[x, 2][[1]], infp[1])
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(pandl[x, 3][[1]])
    crv <- numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      A = matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2,
                 byrow = T)
      b = -c(int[i], int[i + 1])
      crp[i] = solve(A, b)[1]
      crv[i] = solve(A, b)[2]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  laoFin <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(min, pandl[x, 2][[1]], infp[1])
    crp <- c(min, pandl[x, 2][[1]])
    crv <- V(crp)
    tg = int = numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      tg[i] = (V(usepoint[i + 1]) - V(usepoint[i]))/(usepoint[i +
                                                                1] - usepoint[i])
      int[i] = V(usepoint[i]) - tg[i] * usepoint[i]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp[-1]
    result$crv <- crv[-1]
    result
  }
  rtuInf <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(tail(infp, 1), pandl[x, 2][[1]])
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(pandl[x, 3][[1]])
    crv <- numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      A = matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2,
                 byrow = T)
      b = -c(int[i], int[i + 1])
      crp[i] = solve(A, b)[1]
      crv[i] = solve(A, b)[2]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  rtuFin <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(tail(infp, 1), pandl[x, 2][[1]], max)
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(pandl[x, 3][[1]])
    crv <- numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      A = matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2,
                 byrow = T)
      b = -c(int[i], int[i + 1])
      crp[i] = solve(A, b)[1]
      crv[i] = solve(A, b)[2]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  raoFin <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(tail(infp, 1), pandl[x, 2][[1]], max)
    crp <- c(pandl[x, 2][[1]], max)
    crv <- V(crp)
    tg = int = numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      tg[i] = (V(usepoint[i + 1]) - V(usepoint[i]))/(usepoint[i +
                                                                1] - usepoint[i])
      int[i] = V(usepoint[i]) - tg[i] * usepoint[i]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp[-length(crp)]
    result$crv <- crv[-length(crv)]
    result
  }
  ao <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(infp[x - 1], pandl[x, 2][[1]], infp[x])
    crp <- pandl[x, 2][[1]]
    crv <- V(crp)
    tg = int = numeric(pandl[x, 4][[1]])
    for (i in 1:pandl[x, 4][[1]]) {
      tg[i] = (V(usepoint[i + 1]) - V(usepoint[i]))/(usepoint[i +
                                                                1] - usepoint[i])
      int[i] = V(usepoint[i]) - tg[i] * usepoint[i]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  tu <- function(x) {
    result <- list()
    tg = int = crv = crp = c()
    usepoint <- c(infp[x - 1], pandl[x, 2][[1]], infp[x])
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(pandl[x, 3][[1]])
    crv <- numeric(pandl[x, 3][[1]])
    for (i in 1:pandl[x, 3][[1]]) {
      A = matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2,
                 byrow = T)
      b = -c(int[i], int[i + 1])
      crp[i] = solve(A, b)[1]
      crv[i] = solve(A, b)[2]
    }
    result$tg <- tg
    result$int <- int
    result$crp <- crp
    result$crv <- crv
    result
  }
  x_final <- numeric(n)
  for (N in 1:n) {
    p <- function(x) {
      #eval(parse(text = formula))
      fun(x)
    }
    V <- function(x) {
      -log(p(x))
    }
    u = 0
    rate = -1
    while (u >= rate) {
      allpt <- sort(c(sp, infp))
      deriv2 <- function(x) {
        #eval(D(D(parse(text = paste("-log(", formula,
                                   # ")", sp = "")), "x"), "x"))
        d2Minus_Log_Fun(x)
      }
      deriv1 <- function(x) {
       # eval(D(parse(text = paste("-log(", formula, ")",
                                  #sp = "")), "x"))
        dMinus_Log_Fun(x)
      }
      corc <- numeric(length(infp) + 1)
      if (deriv2(infp[length(infp)] + m) > 0 & max == Inf)
        corc[length(infp) + 1] = "rtuInf"
      if (deriv2(infp[length(infp)] + m) > 0 & max != Inf)
        corc[length(infp) + 1] = "rtuFin"
      if (deriv2(infp[length(infp)] + m) < 0 & max != Inf)
        corc[length(infp) + 1] = "raoFin"
      if (deriv2(infp[1] - m) > 0 & min == -Inf)
        corc[1] = "ltuInf"
      if (deriv2(infp[1] - m) > 0 & min != -Inf)
        corc[1] = "ltuFin"
      if (deriv2(infp[1] - m) < 0 & min != -Inf)
        corc[1] = "laoFin"
      if (length(infp) > 1) {
        for (i in 2:length(infp)) {
          if (deriv2(infp[i - 1] + m) < 0)
            corc[i] = "ao"
          if (deriv2(infp[i - 1] + m) > 0)
            corc[i] = "tu"
        }
      }
      parsp <- list()
      parsp[[1]] = sp[which(sp < infp[1])]
      parsp[[length(infp) + 1]] = sp[which(sp > infp[length(infp)])]
      if (length(infp) > 1) {
        for (i in 2:length(infp)) {
          parsp[[i]] = sp[which(sp > infp[i - 1] & sp <
                                  infp[i])]
        }
      }
      pandl <- cbind(corc, parsp, pt = numeric((length(corc))),
                     lns = numeric((length(corc))))
      for (i in 1:nrow(pandl)) {
        if (pandl[i, 1][[1]] == "ao") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]])
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            1
        }
        if (pandl[i, 1][[1]] == "tu") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]]) +
            1
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            2
        }
        if (pandl[i, 1][[1]] == "ltuInf") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]])
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            1
        }
        if (pandl[i, 1][[1]] == "ltuFin") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]]) +
            1
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            2
        }
        if (pandl[i, 1][[1]] == "laoFin") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]]) +
            1
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            1
        }
        if (pandl[i, 1][[1]] == "rtuInf") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]])
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            1
        }
        if (pandl[i, 1][[1]] == "rtuFin") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]]) +
            1
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            2
        }
        if (pandl[i, 1][[1]] == "raoFin") {
          pandl[i, 3][[1]] = length(pandl[i, 2][[1]]) +
            1
          pandl[i, 4][[1]] = length(pandl[i, 2][[1]]) +
            1
        }
      }
      tg_total <- c()
      int_total <- c()
      crp_total <- c()
      crv_total <- c()
      for (i in 1:nrow(pandl)) {
        if (pandl[i, 1][[1]] == "ao")
          tsf <- ao(i)
        if (pandl[i, 1][[1]] == "tu")
          tsf <- tu(i)
        if (pandl[i, 1][[1]] == "raoFin")
          tsf <- raoFin(i)
        if (pandl[i, 1][[1]] == "rtuFin")
          tsf <- rtuFin(i)
        if (pandl[i, 1][[1]] == "rtuInf")
          tsf <- rtuInf(i)
        if (pandl[i, 1][[1]] == "laoFin")
          tsf <- laoFin(i)
        if (pandl[i, 1][[1]] == "ltuFin")
          tsf <- ltuFin(i)
        if (pandl[i, 1][[1]] == "ltuInf")
          tsf <- ltuInf(i)
        tg_total <- c(tg_total, tsf$tg)
        int_total <- c(int_total, tsf$int)
        crp_total <- c(crp_total, tsf$crp)
        crv_total <- c(crv_total, tsf$crv)
      }
      fdtfram <- rbind(cbind(crp_total, crv_total), c(min,
                                                      V(min)), c(max, V(max)))
      fdtfram <- rbind(fdtfram, matrix(c(infp, V(infp)),
                                       nrow = length(infp), byrow = F))
      fdtfram <- fdtfram[order(fdtfram[, 1]), ]
      intsum <- c()
      for (i in 1:length(tg_total)) {
        integ <- function(x) {
          exp(-(tg_total[i] * x + int_total[i]))
        }
        intsum[i] <- integrate(integ, fdtfram[i], fdtfram[i +
                                                            1, 1])[[1]]
      }
      rdm <- runif(1)
      cum = c(0, cumsum(intsum/sum(intsum)))
      idx <- which(rdm < cumsum(intsum/sum(intsum)))[1]
      ifelse(idx > 1, x_star <- (log(-(rdm - cum[idx]) *
                                       sum(intsum) * tg_total[idx] + exp(-tg_total[idx] *
                                                                           fdtfram[idx, 1][[1]] - int_total[idx])) + int_total[idx])/(-tg_total[idx]),
             x_star <- (log(-rdm * sum(intsum) * tg_total[1] +
                              exp(-tg_total[1] * fdtfram[2, 1][[1]] - int_total[1])) +
                          int_total[1])/(-tg_total[1]))
      u <- runif(1)
      rate <- p(x_star)/exp(-tg_total[idx] * x_star - int_total[idx])
      sp = sort(c(sp, x_star))
    }
    x_final[N] <- x_star
  }
  x_final
}
