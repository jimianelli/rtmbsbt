### Explanation

1. **Data and Parameters**: The `data` and `par` lists contain the input data and parameters, respectively. These are passed to the `objective_function`.
2. **Objective Function**: The `objective_function` computes the negative log-likelihood (`nll`) and returns a list of results, including likelihood components, predictions, and residuals.
3. **RTMB Compatibility**: The model is written in R but uses RTMB to handle the C++ backend, making it compatible with TMB's optimization and simulation tools.
4. **Report**: The `report` list contains all the outputs, including likelihood components, predictions, and residuals.

This setup allows you to use the model in R while leveraging the performance of TMB's C++ backend. You can optimize the model using `RTMB::MakeADFun` and `RTMB::nlminb`.
