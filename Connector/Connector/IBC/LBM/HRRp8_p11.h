qOut  = Qinfields[p];
qeq   = Qneqfields[p];
offeq = qOut[indR] - qeq[indR];
a1pr_1 = a1pr_1 + H2H3[(1-1)*nvars+p]*offeq;
a1pr_2 = a1pr_2 + H2H3[(2-1)*nvars+p]*offeq;
a1pr_4 = a1pr_4 + H2H3[(4-1)*nvars+p]*offeq;
a1pr_5 = a1pr_5 + H2H3[(5-1)*nvars+p]*offeq;
a1pr_9 = a1pr_9 + H2H3[(9-1)*nvars+p]*offeq;
