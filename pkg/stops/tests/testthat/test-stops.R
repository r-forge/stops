
library(stops)
data(BankingCrisesDistances)
set.seed(210485)
optstrain <- cops(BankingCrisesDistances[,1:69],loss="strain",verbose=1)
optstress <- cops(BankingCrisesDistances[,1:69],loss="stress",verbose=3)
optsammon <- cops(BankingCrisesDistances[,1:69],loss="sammon",verbose=3)
optelastic <- cops(BankingCrisesDistances[,1:69],loss="elastic",verbose=3)
optsstress <- cops(BankingCrisesDistances[,1:69],loss="sstress",verbose=3)
optrstress <- cops(BankingCrisesDistances[,1:69],loss="rstress",verbose=3)
optpowerstress <- cops(BankingCrisesDistances[,1:69],loss="powerstress",verbose=3)
optpowersammon <- cops(BankingCrisesDistances[,1:69],loss="powersammon",verbose=3)
optpowerelastic <- cops(BankingCrisesDistances[,1:69],loss="powerelastic",verbose=3)

optstrain
optstress
optsammon
optelastic
optsstress
optrstress
optpowerstress
optpowersammon
optpowerelastic



plot(optstrain)
plot(optstress)
plot(optsammon)
plot(optelastic)
plot(optsstress)p
lot(optrstress)
plot(optpowerstress)
plot(optpowersammon)
plot(optpowerelastic)

