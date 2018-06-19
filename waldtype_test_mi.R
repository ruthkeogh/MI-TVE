#between-imputation variance matrix
B=matrix(0,nrow=10,ncol=10)
for(imp in 1:nimp){
  
  temp=all.coef[imp,]-pool.coef
  B.element=temp%*%t(temp)
  B=B+B.element
}
B=B/(nimp-1)

#mean within imputation variance

U=matrix(0,nrow=10,ncol=10)
for(imp in 1:nimp){
  U=U+all.var[imp,,]
}
U=U/nimp

#rm

num.param=4
tr.x1=sum(diag(B[c(3,5,7,9),c(3,5,7,9)]%*%solve(U[c(3,5,7,9),c(3,5,7,9)])))
tr.x2=sum(diag(B[c(4,6,8,10),c(4,6,8,10)]%*%solve(U[c(4,6,8,10),c(4,6,8,10)])))

rm.x1=(1+(1/nimp))*tr.x1/num.param
rm.x2=(1+(1/nimp))*tr.x2/num.param

#test stat

d.x1=t(pool.coef[c(3,5,7,9)])%*%solve(U[c(3,5,7,9),c(3,5,7,9)])%*%
  pool.coef[c(3,5,7,9)]/(num.param*(1+rm.x1))

d.x2=t(pool.coef[c(4,6,8,10)])%*%solve(U[c(4,6,8,10),c(4,6,8,10)])%*%
  pool.coef[c(4,6,8,10)]/(num.param*(1+rm.x2))

v=num.param*(nimp-1)

df.x1=4+(v-4)*((1+(1-2/v)*(1/rm.x1))^2)
df.x2=4+(v-4)*((1+(1-2/v)*(1/rm.x2))^2)

pval.x1=1-pf(d.x1,df1=num.param,df2=df.x1)
pval.x2=1-pf(d.x2,df1=num.param,df2=df.x2)