from math import log, exp


def get_empty_list(count: int) -> list:
    return [0 for _ in range(count)]


def get_empty_matrix(height: int, width: int) -> list:
    return [[0 for j in range(width)] for i in range(height)]


element_name = 'Питан'
k = 28
k += 1
p = 14
# k0 = [0.1586, 0.0912, 0.157, 0.165, 0.0945, 0.2375]
# m = [0, 92, 94, 95, 96, 97, 98, 100]
k0 = [0.0085, 0.00016, 1.5e-9, 0.0035]
k0 = [0, *k0]
k0_amount = sum(k0)
m = [235, 234, 232, 236, 238]
m = [0, *m]
count_m = len(m)
f = 1.3

xi = get_empty_matrix(count_m, 201)
cp = get_empty_matrix(count_m, 201)
cm = get_empty_matrix(count_m, 201)
c = get_empty_matrix(count_m, 201)

ccm = get_empty_matrix(count_m, k)
ccp = get_empty_matrix(count_m, k)
cc = get_empty_matrix(count_m, k)

gp = get_empty_list(k)
gm = get_empty_list(k)
g = get_empty_list(k)
tet = get_empty_list(k)

ggm = get_empty_list(k)
ggp = get_empty_list(k)
gg = get_empty_list(k)
fi = get_empty_list(count_m)
a = get_empty_list(count_m)

for i in range(1, k):
    xi[1][i] = f
    for j in range(2, count_m):
        xi[j][i] = exp((m[5] - m[j]) / (m[5] - m[1]) * log(xi[1][i]))

for i in range(1, count_m - 1):
    print(F'xi{i}=', xi[i][1], F'eps{i}5={log(xi[i][1])}')

print()
f = xi[1][1]
sig = f ** (-0.5)

print(F'sigma={sig}')
print()

for i in range(1, count_m):
    fi[i] = sig * xi[i][1] / (1 + sig * xi[i][1])

fi[-1] = sig / (1 + sig)

for i in range(1, count_m):
    print(f'fi{i}={fi[i]}')

print()

for i in range(1, count_m):
    if abs(fi[i] - 0.5) > 1e-15:
        a[i] = 1 / (2 * fi[i] - 1)
    else:
        a[i] = 1e-15
    print(f'a{i}={a[i]}')

print()

l = [fi[i] / (1 - fi[i]) for i in range(1, count_m)]
l = [0, *l]
for i in range(1, count_m):
    print(f'l{i}={l[i]}')

print()

f = get_empty_list(count_m)
f_amount = 0
for i in range(1, count_m - 1):
    if abs(l[i] - 1) > 1e-10:
        f[i] = k0[i] / (exp(p * log(l[i])) - exp((p - k - 1) * log(l[i]))) * (exp(p * log(l[i])) - 1)
    else:
        f[i] = k0[i] * p / (k + 1)

f[-1] = (1 - k0_amount) / (exp(p * log(l[count_m - 1])) - exp((p - k - 1) * log(l[count_m - 1]))) * (
        exp(p * log(l[count_m - 1])) - 1)
f_amount = sum(f)

for i in range(1, count_m):
    cp[i][k] = f[i] / f_amount

f = get_empty_list(count_m)

for i in range(1, count_m - 1):
    if abs(l[i] - 1) > 1e-10:
        f[i] = k0[i] / (exp(p * log(l[i])) - exp((p - k - 1) * log(l[i]))) * (1 - exp((p - k - 1) * log(l[i])))
    else:
        f[i] = k0[i] * (k - p + 1) / (k + 1)

f[-1] = (1 - k0_amount) / (exp(p * log(l[count_m - 1])) - exp((p - k - 1) * log(l[count_m - 1]))) * (
        1 - exp((p - k - 1) * log(l[count_m - 1])))

f_amount = sum(f)

for i in range(1, count_m):
    cm[i][1] = f[i] / f_amount

kp1 = cp[1][k]
km1 = cm[1][1]
teta = (k0[1] - km1) / (kp1 - km1)

T0 = 100 / 8.76 / 3.6
T = T0 * (k0[1] - km1) / (kp1 - km1)
T1 = T0 - T

print(F'teta={teta} T={T} T1={T1} T0={T0}')

for i in range(1, count_m):
    print(F'cm[{i}, 1]={cm[i][1]}')

print()

for i in range(1, count_m):
    print(F'cp[{i}, k]={cp[i][k]}')

b = get_empty_list(count_m)

for i in range(count_m - 1):
    b[i] = T0 * k0[i] - T * cp[i][k] - T1 * cm[i][1]

cp_amount = 0
cm_amount = 0

for i in range(count_m - 1):
    cp_amount += cp[i][k]
    cm_amount += cm[i][1]

b[-1] = T0 * (1 - k0_amount) - T * (1 - cp_amount) - T1 * (1 - cm_amount)

print('Балансы: ')
for i in range(1, count_m):
    print(f'{i}) {b[i]}')

k_last_f = 1 - k0_amount
k_last_p = 1 - cp_amount
k_last_m = 1 - cm_amount

print(f'Концентрация самого тяжелого компонента:')
print(f'{element_name} {k_last_f}  отб. {k_last_p} отв.  {k_last_m}')

fa = get_empty_list(count_m)

for i in range(1, k):
    for j in range(1, count_m - 1):
        if abs(l[j] - 1) > 1e-10:
            if i >= p:
                fa[j] = T * cp[j][k] * a[j] * (1 - exp((i - k - 1) * log(l[j])))
            else:
                fa[j] = T1 * cm[j][1] * a[j] * (exp(i * log(l[j])) - 1)
        else:
            if i >= p:
                fa[j] = 2 * T * cp[j][k] * (k - i + 1)
            else:
                fa[j] = 2 * T1 * cm[j][1] * i

    if i >= p:
        fa_last = T * (1 - cp_amount) * a[count_m - 1] * (1 - exp((i - k - 1) * log(l[count_m - 1])))
    else:
        fa_last = T1 * (1 - cm_amount) * a[count_m - 1] * (exp(i * log(l[count_m - 1])) - 1)
    fa.append(fa_last)
    fa = [0, *fa]
    f_amount = sum(fa)
    for index in range(1, count_m - 1):
        c[index][i] = fa[index] / f_amount

    fi_fa = get_empty_list(count_m)
    for z in range(1, count_m):
        fi_fa[z] = fa[z] * fi[z]

    f_amount = sum(fi_fa)

    gp[i] = f_amount

    for z in range(1, count_m - 1):
        cp[z][i] = fi_fa[z] / f_amount

    for z in range(1, count_m):
        fi_fa[z] = fa[z] - fi_fa[z]

    f_amount = sum(fi_fa)

    gm[i] = f_amount

    for z in range(1, count_m - 1):
        cm[z][i] = fi_fa[z] / f_amount

    g[i] = gp[i] + gm[i]
    tet[i] = gp[i] / g[i]

for i in range(1, count_m - 1):
    for j in range(1, k):
        ccm[i][j] = cm[i][j]
        ccp[i][j] = cp[i][j]
        cc[i][j] = c[i][j]

for i in range(1, k):
    ggm[i] = gm[i]
    ggp[i] = gp[i]
    gg[i] = g[i]

"""CON"""
for i in range(1, count_m - 1):
    print(F'{element_name}-{m[i]}')
    for j in range(1, k):
        print(F'cm[{i},{j}]={ccm[i][j]} cp[{i},{j}]={ccp[i][j]} c[{i},{j}]={cc[i][j]}')

print(F'{element_name}-{m[-1]}')
for i in range(1, k):
    ccm_amount = 0
    ccp_amount = 0
    cc_amount = 0
    for j in range(1, count_m - 1):
        ccm_amount += ccm[j][i]
        ccp_amount += ccp[j][i]
        cc_amount += cc[j][i]
    print(
        F"cm[{count_m - 1},{i}]={1 - ccm_amount} cp[{count_m - 1}, {i}]={1 - ccp_amount} c[{count_m - 1}, {i}]={1 - cc_amount}")

print('SIGMA')
for i in range(1, k):
    z1, z2, z3, z4 = ggp[i], ccp[1][i], gg[i], cc[1][i]
    f_ = ggp[i] * ccp[1][i] / gg[i] / cc[1][i]
    f1 = f_ / xi[1][i] / (1 - f_)
    f_ = ggp[i] * ccp[2][i] / gg[i] / cc[2][i]
    f2 = f_ / xi[2][i] / (1 - f_)
    print(F"sig[1,{i}]={f1} sig[2,{i}]={f2}")
