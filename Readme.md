# Readme

This project corresponds to the paper **"A Two-Stage Profit Allocation Method for Virtual Power Plants Considering Uncertainties"**. This repository provides partial code and interface structures for our research project. It aims to:

- Show the core implementation structure.
- Help reviewers or interested researchers understand or partially reproduce the results;
- Serve as a basis for full open-sourcing in the future.
- Update for the 10,50 member case. The detailed illustration is in the readme of case4-4

❗ This project is finished and submitted.  Please do not reuse or redistribute without permission. 

##  Project Structure

```
├── README.md  
│   └── Project documentation.

├── Case_4_1_Energy_Market_Toy_Case/                     # Toy case for the energy market  
│   ├── Energy_Market_Case.m                             # Code for optimization and value computation  
│   ├── Allocation.m                             # Allocation_Solution_of the proposed method
│   └── Case4-1-1.xlsx                  # Profit allocation solution of the proposed method  
│   └── Case4-1-2.xlsx                  # Profit allocation process of the proposed allocation method  
│   └── Fig_Three_Member                #Figure Plot for Figure 5:Individual and Allocated Profit for members
├── Case_4_2 Case/                              # Real-world case for the  demand response market  
│   ├── Allocation_Analysis.m                                  # Code for profit creation and calculation  
│   └── Case_Fig_7.xlsx                                # Data for figure plot Figure 7

├── Case_4_3_Case/                       # case (based on actual operation data)
│   └── Parameter/                  #  the detailed parameters of all resources.(Ag:the )
│   └── Allocation_Solution_Paper_Data/                  # Allocation results used in the paper  
│       ├── Without_operator/                            # Data used
│       └── Consider_Operator/                           # Allocation solutions including operator   
│   └── Optimization_Solution/                  # Optimization results used in the paper  
│       ├── Result_DA/                            # Day-Ahead Expectaion
│       ├── Result_RT /                   # Real-Time Operation Solution
│       ├── Result_Mar/                  # Actual Real-Time Operation
│    *Note: Due to the sensitivity of real operational data, full data access is restricted. Please contact |			us if needed.However, the profit allocation results can still be reproduced using the provided   |           matrices in  without_Operator 
│   └── Deviation_Check                  # Check the value of the deviation between the expectation and the |                                         ex-post average
│   └── Case_4_3_Figure_8/                # Plot the figure for figure 8
├── Case_4_4_Case/                       # enlarged case (based on actual operation data)
|
|
|
├── Allocation_Functions/                                # Core functions used for profit allocation  
│   ├── compute_Least_Core.m                             # Least Core allocation implementation (Simplified Nucleolus-First Step)
│   ├── compute_Nucleolus.m                             # Nucleolus allocation implementation (calculate the |                                                           allocation of the nucleolus method)
│   ├── compute_dp.m                                     # Equal DP allocation method  
│   ├── compute_shapley.m                                # Shapley value allocation  
│   ├── Operator_2.m                                     # Matrix generation with operator consideration  
│   ├── Binary_list.m                                    # Generates all member combinations  
│   ├── Consider_Operator.m                              # Allocation calculation considering the operator  
│   └── Without_Operator.m                               # Matrix generation for the "Consider_Operator" function  
```

## Project Status

This project is currently submitted.  
The associated paper is in the writing stage and has not yet been submitted for peer review.  
We plan to make the full code, data, and documentation publicly available after acceptance.

#  Notes

To protect intellectual property, certain parts of the code or data may be redacted or simplified.
 If you require further clarification, feel free to contact us via the email : wsc23@mails.tsinghua.edu.cn

##  License

This code is currently released for academic preview only and is not yet licensed for public use.  
Please do not reuse or redistribute without permission.  
Final licensing terms will be added after formal publication.
