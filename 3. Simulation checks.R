#Checks

prop.table(table(D$family_history_allergyii))
prop.table(table(D$sexii))
prop.table(table(D$mother_asthmaii))
prop.table(table(D$father_asthmaii))
prop.table(table(D$parent_cob_3cat))
prop.table(table(D$number_sibs))
prop.table(table(D$SEIFA_quint))



model <- lm(mother_age_atbirth~as.factor(SEIFA_quint),data=D)
summary(model)

model <- glm(smoke_mother_pregnancyii ~as.factor(SEIFA_quint),data=D,family=binomial(logit))
summary(model)



model <- glm(pre_term ~mother_age_atbirth+as.factor(SEIFA_quint)+smoke_mother_pregnancyii,data=D,family=binomial(logit))
summary(model)



model <- glm(mode_delivery  ~mother_age_atbirth+as.factor(SEIFA_quint),data=D,family=binomial(logit))
summary(model)



model <- glm(LBW_infants ~mode_delivery+pre_term+ mother_age_atbirth+as.factor(SEIFA_quint)+
               smoke_mother_pregnancyii,data=D,family=binomial(logit))
summary(model)



model <- glm(bf_true~sexii+as.factor(number_sibs)+as.factor(parent_cob_3cat)+LBW_infants+mother_asthmaii+
               father_asthmaii+ family_history_allergyii+
               mode_delivery+pre_term+ mother_age_atbirth+as.factor(SEIFA_quint)+
               smoke_mother_pregnancyii,data=D,family=binomial(logit))
summary(model)



model <-  glm(asthma_dx_current_6y~bf_true+S+bf_true*S+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
                LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                smoke_mother_pregnancyii, data=D, family=binomial(logit))

summary(model)

exp(coefficients(model))


#0.85 OR (in data gen 0.80)
#interaction: 0.19  (in data gen 0.2)
