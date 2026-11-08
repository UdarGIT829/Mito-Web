class_name GroupTests
extends RefCounted
## Classical tests for comparing groups of patients (or any units) on one
## measure, matching Stata's defaults so published tables reproduce:
##   kruskal_wallis   Stata kwallis: H with tie correction, chi-squared p
##   rank_sum         Stata ranksum: normal approximation with tie
##                    correction (exact enumeration optional)
##   fisher_exact     r x c exact test (Freeman-Halton), as Stata tab, exact
##   anova            one-way ANOVA F test
##   spearman         rho with Stata's t-based p
## All p-values are two-sided.


## Average ranks (1-based) with ties sharing the mean rank.
static func ranks(values: Array) -> PackedFloat64Array:
	var order := range(values.size())
	order.sort_custom(func(a, b): return values[a] < values[b])
	var out := PackedFloat64Array()
	out.resize(values.size())
	var i := 0
	while i < order.size():
		var j := i
		while j + 1 < order.size() and values[order[j + 1]] == values[order[i]]:
			j += 1
		var mean_rank := (i + j) / 2.0 + 1.0
		for k in range(i, j + 1):
			out[order[k]] = mean_rank
		i = j + 1
	return out


## Sum of t^3 - t over tie groups.
static func _tie_term(values: Array) -> float:
	var counts := {}
	for v in values:
		counts[v] = counts.get(v, 0) + 1
	var t := 0.0
	for c in counts.values():
		t += float(c) * c * c - c
	return t


## Kruskal-Wallis across groups (Arrays of numbers). Returns {"h", "df", "p"}.
static func kruskal_wallis(groups: Array) -> Dictionary:
	var used := groups.filter(func(g): return g.size() > 0)
	var all: Array = []
	for g in used:
		all.append_array(g)
	var n := all.size()
	if used.size() < 2 or n < 3:
		return {"h": NAN, "df": 0, "p": NAN}
	var r := ranks(all)
	var h := 0.0
	var i := 0
	for g in used:
		var s := 0.0
		for k in g.size():
			s += r[i + k]
		i += g.size()
		h += s * s / g.size()
	h = 12.0 / (n * (n + 1.0)) * h - 3.0 * (n + 1.0)
	var correction := 1.0 - _tie_term(all) / (float(n) * n * n - n)
	if correction <= 0.0:
		return {"h": NAN, "df": used.size() - 1, "p": NAN}  # all values equal
	h /= correction
	var df := used.size() - 1
	return {"h": h, "df": df, "p": chi2_sf(h, df)}


## Wilcoxon rank-sum (Mann-Whitney) for two groups. Returns {"z", "p"}.
## exact = true enumerates every split (small groups only).
static func rank_sum(a: Array, b: Array, exact: bool = false) -> Dictionary:
	var n1 := a.size()
	var n2 := b.size()
	if n1 == 0 or n2 == 0:
		return {"z": NAN, "p": NAN}
	var all := a + b
	var n := n1 + n2
	var r := ranks(all)
	var t := 0.0
	for k in n1:
		t += r[k]
	var expected := n1 * (n + 1.0) / 2.0
	if exact and _choose(n, n1) <= 2.0e6:
		return {"z": NAN, "p": _exact_rank_sum(r, n1, absf(t - expected))}
	var variance := n1 * n2 / 12.0 * ((n + 1.0) - _tie_term(all) / (n * (n - 1.0)))
	if variance <= 0.0:
		return {"z": NAN, "p": NAN}
	var z := (t - expected) / sqrt(variance)
	return {"z": z, "p": 2.0 * normal_sf(absf(z))}


static func _exact_rank_sum(r: PackedFloat64Array, n1: int, observed: float) -> float:
	var n := r.size()
	var expected := n1 * (n + 1.0) / 2.0
	var hits := [0, 0]  # [extreme, total]
	_combinations(r, n1, 0, 0.0, expected, observed - 1e-9, hits)
	return float(hits[0]) / hits[1]


static func _combinations(r: PackedFloat64Array, left: int, start: int, sum: float, expected: float, threshold: float, hits: Array) -> void:
	if left == 0:
		hits[1] += 1
		if absf(sum - expected) >= threshold:
			hits[0] += 1
		return
	for i in range(start, r.size() - left + 1):
		_combinations(r, left - 1, i + 1, sum + r[i], expected, threshold, hits)


## Fisher's exact test for an r x c table of counts (Arrays of ints).
## p = sum of probabilities of tables with the same margins that are no
## more likely than the observed one. Returns {"p"}.
static func fisher_exact(table: Array) -> Dictionary:
	var rows: Array = table.filter(func(row): return row.reduce(func(s, x): return s + x, 0) > 0)
	if rows.size() < 2:
		return {"p": NAN}
	var cols: int = rows[0].size()
	var col_sums: Array = []
	for c in cols:
		col_sums.append(rows.reduce(func(s, row): return s + row[c], 0))
	var keep: Array = range(cols).filter(func(c): return col_sums[c] > 0)
	if keep.size() < 2:
		return {"p": 1.0}
	var trimmed: Array = rows.map(func(row): return keep.map(func(c): return row[c]))
	var row_sums: Array = trimmed.map(func(row): return row.reduce(func(s, x): return s + x, 0))
	col_sums = keep.map(func(c): return col_sums[c])
	var n: int = row_sums.reduce(func(s, x): return s + x, 0)
	var log_const := -_log_factorial(n)
	for s in row_sums:
		log_const += _log_factorial(s)
	for s in col_sums:
		log_const += _log_factorial(s)
	var observed := log_const
	for row in trimmed:
		for x in row:
			observed -= _log_factorial(x)
	# Exact enumeration when the number of tables is manageable; otherwise a
	# seeded Monte Carlo over random tables with the same margins.
	var tables := 1.0
	for i in trimmed.size() - 1:
		tables *= _choose(row_sums[i] + col_sums.size() - 1, col_sums.size() - 1)
	if tables <= MAX_EXACT_TABLES:
		var acc := [0.0]
		_fisher_walk(row_sums, col_sums.duplicate(), 0, log_const, observed + 1e-7, acc)
		return {"p": minf(1.0, acc[0]), "method": "exact"}
	return {"p": _fisher_monte_carlo(row_sums, col_sums, log_const, observed + 1e-7), "method": "Monte Carlo"}


const MAX_EXACT_TABLES := 3.0e6
const MONTE_CARLO_TABLES := 20000


## Random tables with fixed margins: shuffle the column labels of the units.
static func _fisher_monte_carlo(row_sums: Array, col_sums: Array, log_const: float, limit: float) -> float:
	var rng := RandomNumberGenerator.new()
	rng.seed = 20260101
	var labels := PackedInt32Array()
	for c in col_sums.size():
		for _k in col_sums[c]:
			labels.append(c)
	var hits := 0
	for _t in MONTE_CARLO_TABLES:
		for i in range(labels.size() - 1, 0, -1):
			var j := rng.randi_range(0, i)
			var tmp := labels[i]
			labels[i] = labels[j]
			labels[j] = tmp
		var lp := log_const
		var at := 0
		for r in row_sums:
			var counts := {}
			for _k in r:
				counts[labels[at]] = counts.get(labels[at], 0) + 1
				at += 1
			for x in counts.values():
				lp -= _log_factorial(x)
		if lp <= limit:
			hits += 1
	return (hits + 1.0) / (MONTE_CARLO_TABLES + 1.0)


## Enumerates tables row by row: each row splits its total across the
## remaining column capacities.
static func _fisher_walk(row_sums: Array, cols_left: Array, row: int, log_p: float, limit: float, acc: Array) -> void:
	if row == row_sums.size() - 1:
		var lp := log_p
		for c in cols_left:
			lp -= _log_factorial(c)
		if lp <= limit:
			acc[0] += exp(lp)
		return
	_fisher_row(row_sums, cols_left, row, 0, row_sums[row], log_p, limit, acc)


static func _fisher_row(row_sums: Array, cols_left: Array, row: int, col: int, remaining: int, log_p: float, limit: float, acc: Array) -> void:
	if col == cols_left.size() - 1:
		if remaining > cols_left[col]:
			return
		cols_left[col] -= remaining
		_fisher_walk(row_sums, cols_left, row + 1, log_p - _log_factorial(remaining), limit, acc)
		cols_left[col] += remaining
		return
	var later := 0
	for c in range(col + 1, cols_left.size()):
		later += cols_left[c]
	for x in range(maxi(0, remaining - later), mini(remaining, cols_left[col]) + 1):
		cols_left[col] -= x
		_fisher_row(row_sums, cols_left, row, col + 1, remaining - x, log_p - _log_factorial(x), limit, acc)
		cols_left[col] += x


## One-way ANOVA. Returns {"f", "df1", "df2", "p"}.
static func anova(groups: Array) -> Dictionary:
	var used := groups.filter(func(g): return g.size() > 0)
	var all: Array = []
	for g in used:
		all.append_array(g)
	var n := all.size()
	var k := used.size()
	if k < 2 or n <= k:
		return {"f": NAN, "df1": k - 1, "df2": n - k, "p": NAN}
	var grand := _mean(all)
	var between := 0.0
	var within := 0.0
	for g in used:
		var m := _mean(g)
		between += g.size() * (m - grand) * (m - grand)
		for x in g:
			within += (x - m) * (x - m)
	if within <= 0.0:
		return {"f": NAN, "df1": k - 1, "df2": n - k, "p": NAN}
	var f := (between / (k - 1)) / (within / (n - k))
	return {"f": f, "df1": k - 1, "df2": n - k, "p": f_sf(f, k - 1, n - k)}


## Spearman rank correlation of paired values. Returns {"rho", "n", "p"}.
static func spearman(x: Array, y: Array) -> Dictionary:
	var n := x.size()
	if n < 3:
		return {"rho": NAN, "n": n, "p": NAN}
	var rx := ranks(x)
	var ry := ranks(y)
	var mx := _mean(Array(rx))
	var my := _mean(Array(ry))
	var sxy := 0.0
	var sxx := 0.0
	var syy := 0.0
	for i in n:
		sxy += (rx[i] - mx) * (ry[i] - my)
		sxx += (rx[i] - mx) * (rx[i] - mx)
		syy += (ry[i] - my) * (ry[i] - my)
	if sxx <= 0.0 or syy <= 0.0:
		return {"rho": NAN, "n": n, "p": NAN}
	var rho := sxy / sqrt(sxx * syy)
	if absf(rho) >= 1.0:
		return {"rho": rho, "n": n, "p": 0.0}
	var t := rho * sqrt((n - 2.0) / (1.0 - rho * rho))
	return {"rho": rho, "n": n, "p": t_sf_two_sided(t, n - 2)}


# --- Summaries ---------------------------------------------------------------

static func _mean(values: Array) -> float:
	var s := 0.0
	for v in values:
		s += v
	return s / values.size() if values.size() else NAN


static func mean_sd(values: Array) -> Array:
	var n := values.size()
	if n == 0:
		return [NAN, NAN]
	var m := _mean(values)
	if n < 2:
		return [m, NAN]
	var ss := 0.0
	for v in values:
		ss += (v - m) * (v - m)
	return [m, sqrt(ss / (n - 1))]


## Median and 25th/75th percentiles, as Stata's centile/summarize detail
## (empirical distribution, averaging at exact boundaries).
static func median_iqr(values: Array) -> Array:
	if values.is_empty():
		return [NAN, NAN, NAN]
	var s := values.duplicate()
	s.sort()
	return [_percentile(s, 50.0), _percentile(s, 25.0), _percentile(s, 75.0)]


static func _percentile(sorted: Array, p: float) -> float:
	var n := sorted.size()
	var pos := n * p / 100.0
	var i := int(floor(pos))
	if absf(pos - i) < 1e-9:
		if i <= 0:
			return float(sorted[0])
		if i >= n:
			return float(sorted[n - 1])
		return (sorted[i - 1] + sorted[i]) / 2.0
	return float(sorted[mini(i, n - 1)])


# --- Distributions -----------------------------------------------------------

static func normal_sf(z: float) -> float:
	return 0.5 * _erfc(z / sqrt(2.0))


static func _erfc(x: float) -> float:
	# Numerical Recipes erfcc: fractional error < 1.2e-7.
	var z := absf(x)
	var t := 1.0 / (1.0 + 0.5 * z)
	var r := t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806
		+ t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
	return r if x >= 0.0 else 2.0 - r


static func chi2_sf(x: float, df: int) -> float:
	if x <= 0.0:
		return 1.0
	return 1.0 - _gamma_p(df / 2.0, x / 2.0)


static func f_sf(f: float, d1: int, d2: int) -> float:
	if f <= 0.0:
		return 1.0
	return _beta_inc(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * f))


static func t_sf_two_sided(t: float, df: int) -> float:
	return _beta_inc(df / 2.0, 0.5, df / (df + t * t))


static func _log_gamma(x: float) -> float:
	var c := [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
	var y := x
	var tmp := x + 5.5
	tmp -= (x + 0.5) * log(tmp)
	var ser := 1.000000000190015
	for coefficient in c:
		y += 1.0
		ser += coefficient / y
	return -tmp + log(2.5066282746310005 * ser / x)


static func _log_factorial(n: int) -> float:
	return _log_gamma(n + 1.0)


static func _choose(n: int, k: int) -> float:
	return exp(_log_factorial(n) - _log_factorial(k) - _log_factorial(n - k))


## Regularized lower incomplete gamma P(a, x).
static func _gamma_p(a: float, x: float) -> float:
	if x < a + 1.0:
		var ap := a
		var del := 1.0 / a
		var sum := del
		for _i in 500:
			ap += 1.0
			del *= x / ap
			sum += del
			if absf(del) < absf(sum) * 1e-14:
				break
		return sum * exp(-x + a * log(x) - _log_gamma(a))
	var b := x + 1.0 - a
	var c := 1.0 / 1e-300
	var d := 1.0 / b
	var h := d
	for i in range(1, 500):
		var an := -i * (i - a)
		b += 2.0
		d = an * d + b
		if absf(d) < 1e-300:
			d = 1e-300
		c = b + an / c
		if absf(c) < 1e-300:
			c = 1e-300
		d = 1.0 / d
		var delta := d * c
		h *= delta
		if absf(delta - 1.0) < 1e-14:
			break
	return 1.0 - exp(-x + a * log(x) - _log_gamma(a)) * h


## Regularized incomplete beta I_x(a, b).
static func _beta_inc(a: float, b: float, x: float) -> float:
	if x <= 0.0:
		return 0.0
	if x >= 1.0:
		return 1.0
	var front := exp(_log_gamma(a + b) - _log_gamma(a) - _log_gamma(b) + a * log(x) + b * log(1.0 - x))
	if x < (a + 1.0) / (a + b + 2.0):
		return front * _beta_cf(a, b, x) / a
	return 1.0 - front * _beta_cf(b, a, 1.0 - x) / b


static func _beta_cf(a: float, b: float, x: float) -> float:
	var qab := a + b
	var qap := a + 1.0
	var qam := a - 1.0
	var c := 1.0
	var d := 1.0 - qab * x / qap
	if absf(d) < 1e-300:
		d = 1e-300
	d = 1.0 / d
	var h := d
	for m in range(1, 300):
		var m2 := 2.0 * m
		var aa := m * (b - m) * x / ((qam + m2) * (a + m2))
		d = 1.0 + aa * d
		if absf(d) < 1e-300:
			d = 1e-300
		c = 1.0 + aa / c
		if absf(c) < 1e-300:
			c = 1e-300
		d = 1.0 / d
		h *= d * c
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
		d = 1.0 + aa * d
		if absf(d) < 1e-300:
			d = 1e-300
		c = 1.0 + aa / c
		if absf(c) < 1e-300:
			c = 1e-300
		d = 1.0 / d
		var delta := d * c
		h *= delta
		if absf(delta - 1.0) < 1e-14:
			break
	return h
