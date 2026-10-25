# Publish to a new standalone GitHub repository

Create a **new** repository named `repair-filtered-mutational-histories` under your GitHub account. It should **not** be placed inside or pushed to your earlier Haldane papers repository.

## GitHub website method

1. GitHub → **New repository** → name: `repair-filtered-mutational-histories`.
2. Choose **Public** once ready for open review. Avoid initializing with another README; this package already includes one.
3. Extract the delivered ZIP into a local directory. Use GitHub Desktop or `git` to commit that directory to the new repository. Do not upload the ZIP as a single file instead of extracting it.
4. Open the README on GitHub and verify that the figure links, folders, and citation metadata display as intended.
5. Add an author-approved open-source `LICENSE` and complete any required disclosure/attribution before public release.

## Terminal / GitHub CLI method (optional)

```bash
cd repair-filtered-mutational-histories
git init -b main
git add .
git commit -m "Release reproducibility code for Paper I"
gh auth login                          # once, if not already authenticated
gh repo create repair-filtered-mutational-histories --public --source=. --remote=origin --push
```

If your environment lacks GitHub CLI, create the empty repo on github.com and then use the remote URL shown by GitHub:

```bash
git remote add origin https://github.com/YOUR-USERNAME/repair-filtered-mutational-histories.git
git push -u origin main
```

The `YOUR-USERNAME` placeholder must be replaced with your actual GitHub username; no URL has been invented or published by this package.

## After publishing

- Run the GitHub Actions reproducibility workflow; investigate any failure before claiming public reproducibility.
- Consider GitHub **Releases** and Zenodo archival to obtain a versioned DOI for the exact scripts accompanying the submitted article.
- Replace the manuscript's data/code availability placeholder with the **real**, publicly accessible GitHub URL (and archived DOI if applicable).
- If later updating the repository, tag the release matching the accepted manuscript so readers can reproduce the exact submitted figures.
