suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(lubridate)
  library(ggplot2)
})

dir.create("figs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 先把两个 8 天的原始文件读进来，按日期合成一张表
era <- fread("ERA5L_TWSA.csv")
gpp <- fread("GPP_anomaly.csv")

era[, date := ymd(date)]
gpp[, date := ymd(date)]

dat <- merge(era[, .(date, VPD, SM1)], gpp[, .(date, GPP)], by = "date")
setorder(dat, date)

cat("merged rows:", nrow(dat), "from", format(min(dat$date)),
    "to", format(max(dat$date)), "\n")

# 标准化成干旱指数（按排名转成正态得分），不然 VPD 和 SM1 单位不一样没法比；
# 这里约定越大越干，SM1 越大其实越湿，所以给它加个负号
standardize <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(x))
}

dat[, drought_met := standardize(VPD)]
dat[, drought_ag := -standardize(SM1)]

cat("cor(drought_met, drought_ag) =", round(cor(dat$drought_met, dat$drought_ag), 3), "\n")

# 连续超阈值的时段就算一次干旱事件，顺手把起点、持续、强度、烈度都算出来
extract_events <- function(data, col, thr) {
  e <- data.table(date = data$date, idx = data[[col]])
  e[, dry := idx > thr]
  e[, id := rleid(dry)]

  ev <- e[dry == TRUE, .(
    start = min(date),
    end = max(date),
    duration = .N,
    intensity = mean(idx),
    severity = sum(idx - thr)
  ), by = id]

  setorder(ev, start)
  ev[, id := NULL]
  ev
}

met_events <- extract_events(dat, "drought_met", 1)
ag_events <- extract_events(dat, "drought_ag", 1)

cat("\nmeteorological drought events:", nrow(met_events),
    " mean duration:", round(mean(met_events$duration), 1), "\n")
cat("agricultural  drought events:", nrow(ag_events),
    " mean duration:", round(mean(ag_events$duration), 1), "\n")

fwrite(met_events, "results/met_events.csv")
fwrite(ag_events, "results/ag_events.csv")

# 用半年窗口做滑动相关，24 期 x 8 天
window <- 24
roll_cor <- rep(NA_real_, nrow(dat))
for (end in seq(window, nrow(dat))) {
  roll_cor[end] <- cor(dat$drought_met[(end - window + 1):end], dat$drought_ag[(end - window + 1):end])
}
dat[, rolling_cor := roll_cor]

# 相关系数高就是“同步”，否则不同步；前 23 天还没攒够窗口是 NA，这里也留 NA
dat[, regime := fifelse(rolling_cor >= 0.5, "synchronous", "asynchronous")]

regime_share <- dat[!is.na(regime), .(n = .N), by = regime]
regime_share[, share := round(n / sum(n), 3)]
cat("\nsynchronous / asynchronous share:\n")
print(regime_share)

# 算 GPP 对干旱的响应，整体和分时段用同一套算法，省得重复抄
gpp_stats <- function(d) d %>% summarise(
  n = n(),
  r_ag = round(cor(GPP, drought_ag), 3),
  r_met = round(cor(GPP, drought_met), 3),
  slope = round(coef(lm(GPP ~ drought_ag))[2], 4),
  dry = round(mean(GPP[drought_ag > 1]), 3),
  wet = round(mean(GPP[drought_ag <= 1]), 3)
)

overall_stats <- dat %>% gpp_stats() %>% mutate(regime = "all")
regime_stats <- dat %>% filter(!is.na(regime)) %>% group_by(regime) %>% gpp_stats() %>% ungroup()

gpp_response <- bind_rows(overall_stats, regime_stats)
cat("\nGPP response by regime:\n")
print(gpp_response)
fwrite(gpp_response, "results/gpp_response.csv")

# 出图
theme_set(theme_bw(base_size = 12))

p1 <- dat[, .(date, drought_met, drought_ag)] %>%
  melt(id.vars = "date", variable.name = "name", value.name = "value") %>%
  ggplot(aes(date, value, color = name)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(drought_met = "#c0392b", drought_ag = "#2471a3"),
                     labels = c(drought_met = "Meteorological (VPD)", drought_ag = "Agricultural (SM)")) +
  labs(x = NULL, y = "Standardized drought index", color = NULL,
       title = "Standardized drought indices (positive = drier)") +
  theme(legend.position = "top")
ggsave("figs/drought_indices.png", p1, width = 9, height = 4, dpi = 150)

p2 <- dat[!is.na(rolling_cor)] %>%
  ggplot(aes(date, rolling_cor)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_line(color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = regime), size = 1) +
  scale_color_manual(values = c(synchronous = "#1e8449", asynchronous = "#b9770e")) +
  labs(x = NULL, y = "24-window rolling correlation", color = NULL,
       title = "Sliding correlation: synchronous vs asynchronous coupling") +
  theme(legend.position = "top")
ggsave("figs/rolling_correlation.png", p2, width = 9, height = 4, dpi = 150)

p3 <- dat[!is.na(regime)] %>%
  ggplot(aes(drought_ag, GPP, color = regime)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ regime) +
  scale_color_manual(values = c(synchronous = "#1e8449", asynchronous = "#b9770e")) +
  labs(x = "Agricultural drought index", y = "GPP anomaly",
       title = "GPP response to agricultural drought by coupling regime") +
  theme(legend.position = "none")
ggsave("figs/gpp_response.png", p3, width = 9, height = 4, dpi = 150)

cat("\nfigures written to figs/ and tables to results/\n")
