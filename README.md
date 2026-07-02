# Production Template

This template should be used when developing an app or pipeline intended for production.

This template enforces restrictions on users who are not the author of the repository to ensure code integrity and security.

### Repository Rules

1. **Require a Pull Request Before Merging:**  
   All commits must be made to a non-protected branch. The `main` branch is protected, and changes must be submitted via a pull request (PR) before merging.

2. **Require Approvals:**  
   Pull requests targeting the `main` branch require at least **one approval** and **no requested changes** before they can be merged.

3. **Require Review from Code Owners:**  
   An approved review is mandatory for pull requests that include files with a designated code owner.

4. **Admin Bypass for Emergencies:**  
   Administrators are allowed to bypass required pull requests in case of emergencies.

5. **Require Signed Commits:**  
   Commits pushed to the `main` branch must have verified signatures to ensure authenticity.

6. **Test Validation via GitHub Actions:**  
   A GitHub Actions workflow is set up to run tests automatically when commits or pull requests are made to the `main` branch. Commits that do not pass the test suite will not be merged.

### Data Handling Policy

This template includes a README template that explicitly prohibits the upload of any genomics or other personally identifiable data to the online repository. Please ensure all data shared in this repository complies with data protection and privacy regulations.
