# This file is part of a program that make prediction of active
# protein phosphorylation sites using machine learning
# Copyright (C) 2018  Zoé Brunet
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import h2o
import matplotlib
import matplotlib.pyplot as plt

h2o.init(nthreads=-1)
mod = h2o.load_model("path")
thresholds = [x/100.0 for x in range(1, 100)]
precision = mod.precision(valid=True, thresholds=thresholds)
recall = mod.metric("recall", valid=True, thresholds=thresholds)
for l in [precision, recall]:
    for i, elem in enumerate(l):
        l[i] = elem[1]
plt.subplot(121)
plt.plot(recall, precision)
plt.xlabel('recall')
plt.ylabel('precision')
plt.title('Precision-Recall Curve')
plt.subplot(122)
plt.plot(thresholds, thresholds)
perf = mod.model_performance(valid=True)
print(mod.confusion_matrix(valid=True))
perf.plot()
