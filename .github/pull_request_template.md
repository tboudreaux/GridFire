<!--
Thank you for contributing to GridFire! 🚀

To help the reviewers understand your changes, please fill out the sections below.
Feel free to delete any sections that are not relevant to your pull request.
-->
# Your PR Title Here

## 🎯 Type of Change

<!-- Please check the box that best describes your PR. -->

- [ ] 🐞 Bug Fix (a non-breaking change that fixes an issue)

- [ ] ✨ New Feature (a non-breaking change that adds functionality)

- [ ] 💥 Breaking Change (a fix or feature that would cause existing functionality to not work as expected)

- [ ] ⚡ Performance Optimization (a change that improves the speed or memory usage of the code)

- [ ] ♻️ Refactor (a code change that neither fixes a bug nor adds a feature)

- [ ] 📚 Documentation Update

- [ ] 🧪 Test Suite Addition/Update

## 🔗 Related Issue

<!-- If this PR addresses a specific issue, please link to it here. -->

<!-- E.g., "Closes #42" or "Resolves #15" -->

- Related to Issue: #

## 📝 Description of Change

<!--
Provide a clear and concise description of the changes.
    Explain the "what" and the "why" of this pull request.
    What was the problem, and how does this PR solve it?
-->


## 🏛️ Architectural & Implementation Details

<!--
Describe the technical details of your implementation.
    What design patterns did you use (e.g., Strategy Pattern for Screening Models)?
    How does this change interact with other parts of the code (e.g., the Engine, Solver, AD taping)?
    Are there any new classes, and what are their responsibilities?
-->


## 🔬 Physics Impact

<!--
Does this PR introduce new physics (e.g., electron screening, reverse reactions)?
    Does it change the physical assumptions of an existing model?
    How does this affect the physical results?
-->


## ✅ Validation & Testing

<!--
Describe how you have validated your changes.
    A PR will not be merged without validation.
    How do the results compare to a known benchmark (e.g., approx8, MESA)?
    Please include tables or plots showing the comparison.
    If this is a bug fix, describe the test case that now passes.
-->

##### Test Configuration:
##### Results:

<!--Example table format for results-->
| Model   | Final H-1 Mass Fraction | Notes                  | 
|---------|-------------------------|------------------------|
| approx8 | 0.1234                  | Matches expected value |


<!-- You can embed images or plots from your validation runs here. -->
## ⏱️ Performance Impact

<!--
If relevant, describe the performance impact of your changes.
    Provide profiling data before and after the change.
    Include metrics like total execution time and number of solver steps.
-->


## ☑️ Pre-Merge Checklist

<!-- Go over all the following points, and put an x in all the boxes that apply. -->

- [ ] My code follows the style guidelines of this project.

- [ ] I have performed a self-review of my own code.

- [ ] I have commented my code, particularly in hard-to-understand areas.

- [ ] I have made corresponding changes to the documentation.

- [ ] My changes generate no new warnings.

- [ ] I have added tests that prove my fix is effective or that my feature works.

- [ ] New and existing unit tests pass locally with my changes.

- [ ] I have confirmed that this PR does not break the AD tape recording process.