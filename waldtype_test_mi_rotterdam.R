#between-imputation variance matrix
B=matrix(0,nrow=num.param,ncol=num.param)
for(imp in 1:nimp){
  
  temp=all.coef[imp,]-pool.coef
  B.element=temp%*%t(temp)
  B=B+B.element
}
B=B/(nimp-1)

#mean within imputation variance

U=matrix(0,nrow=num.param,ncol=num.param)
for(imp in 1:nimp){
  U=U+all.var[imp,,]
}
U=U/nimp

#rm

tr.x=sum(diag(B%*%solve(U)))

rm.x=(1+(1/nimp))*tr.x/num.param

#test stat

d.x=t(pool.coef)%*%solve(U)%*%
  pool.coef/(num.param*(1+rm.x))

v=num.param*(nimp-1)

df.x=4+(v-4)*((1+(1-2/v)*(1/rm.x))^2)

pval=1-pf(d.x,df1=num.param,df2=df.x)
