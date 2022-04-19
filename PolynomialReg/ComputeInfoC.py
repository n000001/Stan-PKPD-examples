
import numpy as np


def waic(log_lik):
    lpd=np.sum(np.log(np.mean(np.exp(log_lik), axis=0)))
    p_waic = np.sum(np.var(log_lik, axis=0, dtype=np.float128))
    elpd_waic=lpd - p_waic
    waic = -2 * elpd_waic
    return dict(p_waic=p_waic, elpd_waic=elpd_waic, waic=waic)

def aic(log_lik, nparams):
    k = nparams
    return 2 * k - 2 * np.max(log_lik)
    
def aicc(log_lik, nparams, ndata):
    k = nparams
    n = ndata
    return aic(log_lik, nparams) + 2 * k * (k + 1) / (n - k - 1)
    
def bic(log_lik, nparams, ndata):
    k = nparams
    n = ndata
    return -2 * np.max(log_lik) + k * np.log(n)
    
def dic(log_lik, uparams):
    pDIC = 2 * (lik_at_mean(log_lik, uparams) - np.mean(log_lik, axis=0))
    DIC = -2 * (lik_at_mean(log_lik, uparams)-pDIC) 
    return DIC

