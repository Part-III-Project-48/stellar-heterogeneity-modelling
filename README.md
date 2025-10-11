# stellar-heterogeneity-modelling

Cambridge Part III Project 48

PHOENIX homepage : [link](https://phoenix.astro.physik.uni-goettingen.de/)
PHOENIX reference paper : [link](https://arxiv.org/abs/1303.5632v2)

NOTE: the ftp links on the PHOENIX website are broken. Converting the format to start with https:// works; see phoenix-grid-creator/README.md for an example
 
# Conventions

## Commit Types

- **`feat`**: A new feature  
- **`fix`**: A bug fix  
- **`docs`**: Documentation-only changes
- **`notes`**: Adding dev-notes / diary entries
- **`style`**: Code style changes (whitespace, formatting, etc. — no code behavior change)  
- **`refactor`**: A code change that neither fixes a bug nor adds a feature  
- **`perf`**: A code change that improves performance  
- **`test`**: Adding or correcting tests  
- **`build`**: Changes to the build system or dependencies (e.g. npm, Makefile)  
- **`ci`**: Changes to CI configuration or scripts (e.g. GitHub Actions, Travis)
- **`chore`**: General maintenance tasks not related to features, fixes, docs, or build

## Branch Names

Branch names can follow the same names, but are formatted like
- **`feat/adding-x-from-y`**
- **`fix-enemymeshes/removing-extraneous-faces`**

(as **`:`**, **`()`** etc would have to be escaped and branch names cannot have whitespace.)

## Example Messages

```bash
feat(bvh): implement initial bounding volume hierarchy generation
fix(plot): correct axis scaling in star visualisation
docs: update README with usage instructions
build: update python_requirements.txt
chore: remove .vscode folder from git tracking (cache)
```

## Sources

This project follows the [Conventional Commits specification v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/#summary).

Some of the commit types listed above are from the [Angular Commit Message Guidelines](https://github.com/angular/angular/blob/22b96b9/CONTRIBUTING.md#-commit-message-guidelines). That link also gives some general guidelines for message formatting that I think are nice.
