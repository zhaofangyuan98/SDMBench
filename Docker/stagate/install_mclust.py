import rpy2.robjects.packages as r
utils = r.importr('utils')
utils.chooseCRANmirror(ind=17)
package_name = 'mclust'
utils.install_packages(package_name)
