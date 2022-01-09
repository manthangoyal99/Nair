
# Define the function for obtaining A/At (it is a parabola)
def fun_A_over_At_parabolic(xi,xf,Am,At):
    return lambda x: 1 + (Am/At-1)*pow(2*(x-xi)/(xf-xi)-1,2.)
