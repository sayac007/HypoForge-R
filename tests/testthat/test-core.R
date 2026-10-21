test_that('identity and lag semantics', {
  x <- data.frame(a=1:10)
  expect_equal(hf_evaluate(hf_spec('a','identity','a'),x),1:10)
  expect_equal(hf_evaluate(hf_spec('l','lag','a',list(periods=2)),x)[3],1)
})
test_that('no overlapping rolling splits', {
  z<-hf_rolling_splits(120,50,10,gap=2)
  expect_true(all(vapply(z,function(s) max(s$train)+2<min(s$test),logical(1))))
})
test_that('known signal improves rolling RMSE', {
  set.seed(14);n<-320;x<-rnorm(n);signal<-rnorm(n);noise<-rnorm(n)
  d<-data.frame(t=seq_len(n),x=x,signal=signal,noise=noise,y=2*x+3*signal+rnorm(n,sd=.5))
  cands<-list(hf_spec('signal','identity','signal'),hf_spec('noise','identity','noise'))
  out<-hf_tournament(d,'y','t','x',cands,c(x='known_future',signal='known_future',noise='known_future'),
                     min_train=100,horizon=20)
  expect_equal(out$summary$feature[1],'signal')
  expect_gt(out$summary$mean_improvement[1],1)
})
