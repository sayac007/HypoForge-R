test_that('causal geometric adstock', {
  expect_equal(hf_geometric_adstock(c(1,0,0,0),alpha=.5,l_max=3,normalize=FALSE),c(1,.5,.25,.125))
  first <- hf_geometric_adstock(c(1,0,0),alpha=.5,l_max=3,normalize=FALSE)
  altered <- hf_geometric_adstock(c(1,0,999),alpha=.5,l_max=3,normalize=FALSE)
  expect_equal(first[1:2],altered[1:2])
})
test_that('Easter and saturation', {
  expect_equal(as.character(hf_easter(2026)),'2026-04-05')
  expect_equal(hf_michaelis_menten(2,half_saturation=2),.5)
  expect_equal(hf_reach_frequency(100,2,half_saturation=2),50)
})
test_that('calendar events', {
  f <- hf_calendar_features(as.Date(c('2026-04-05','2026-03-29')), offsets_weeks=c(-1,0))
  expect_equal(f$`event_easter_w+0`[1],1)
  expect_equal(f$`event_easter_w-1`[2],1)
})
