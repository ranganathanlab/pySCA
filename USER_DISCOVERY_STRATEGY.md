# User Discovery Strategy for pySCA 7.0

## The Problem

If you make the old `pySCA` repository private, users with:
- Bookmarks to `github.com/ranganathanlab/pySCA`
- Old clone URLs
- Documentation links
- Package dependencies pointing to the old repo

**Will NOT automatically find the new version.**

## Solutions

### Option 1: Keep Old Repository Public but Archived (RECOMMENDED)

**Best for user discovery and continuity**

1. **Keep old `pySCA` public** (don't make it private)
2. **Archive the repository** (GitHub: Settings → Archive this repository)
3. **Add prominent notice** in README pointing to new version
4. **Create new `pySCA` repository** for v7.0 (or rename pySCA-dev)

**What users see:**
- Old repo: Shows "This repository has been archived" banner
- README at top: "⚠️ This is pySCA v6.x (archived). For the latest version, see [pySCA v7.0](https://github.com/ranganathanlab/pySCA)"
- Old links still work, but redirect users to new version
- GitHub search still finds both, but new one is active

**Benefits:**
- ✅ Old links/bookmarks still work
- ✅ Users naturally discover new version
- ✅ Clear migration path
- ✅ No broken links
- ✅ Search engines find both versions

### Option 2: Repository Transfer with Redirect

**If you want to use the same repository name:**

1. **Transfer old `pySCA` to a new name** (e.g., `pySCA-legacy`)
2. **Rename `pySCA-dev` to `pySCA`**
3. GitHub automatically redirects old URLs for a period

**What users see:**
- Old URLs redirect to new location (temporary)
- Eventually old URLs may break
- Users need to update their links

**Limitations:**
- ⚠️ Redirects are temporary (GitHub doesn't guarantee permanent redirects)
- ⚠️ Some tools may not follow redirects
- ⚠️ Less clear than archived repo

### Option 3: Keep Both Public, Add Notices

**Maximum visibility:**

1. **Keep old `pySCA` public** (v6.x)
2. **Make `pySCA-dev` public and rename to `pySCA-v7` or keep as `pySCA`**
3. **Add notices in both repositories**

**Old pySCA README:**
```markdown
# pySCA (v6.x - Legacy)

⚠️ **This is the legacy version of pySCA (v6.x)**

**For the latest version, please use [pySCA v7.0](https://github.com/ranganathanlab/pySCA)**

This repository is maintained for historical reference only.
```

**New pySCA README:**
```markdown
# pySCA v7.0

✅ **This is the current version of pySCA**

For the legacy v6.x version, see [pySCA v6.x](https://github.com/ranganathanlab/pySCA-legacy)
```

**Benefits:**
- ✅ Both versions accessible
- ✅ Clear version distinction
- ✅ Users can choose which version
- ✅ No broken links

**Drawbacks:**
- ⚠️ May confuse users (which one to use?)
- ⚠️ Maintenance overhead (two public repos)

## Recommended Approach: Archive Old + New Public Repo

### Step-by-Step:

1. **Update old `pySCA` README** with notice:
   ```markdown
   # pySCA (v6.x - Archived)
   
   ⚠️ **This repository has been archived.**
   
   **For the latest version (v7.0), please visit:**
   [https://github.com/ranganathanlab/pySCA](https://github.com/ranganathanlab/pySCA)
   
   This version (v6.x) is maintained for historical reference only.
   ```

2. **Archive the old repository:**
   - Go to old `pySCA` → Settings → Danger Zone
   - Click "Archive this repository"
   - Repository becomes read-only with archive banner

3. **Rename `pySCA-dev` to `pySCA`** (or create new public `pySCA`)
   - This becomes the active, public repository
   - Users naturally find it via search
   - Old repo links show archive notice + redirect

4. **Add to new `pySCA` README:**
   ```markdown
   # pySCA v7.0
   
   Python 3 implementation of Statistical Coupling Analysis (SCA).
   
   > **Note:** If you're looking for pySCA v6.x, see the [archived repository](https://github.com/ranganathanlab/pySCA-legacy).
   ```

## How Users Will Discover v7.0

### Automatic Discovery:
1. **GitHub Search**: Searching "pySCA" shows both, but new one is active
2. **Repository Links**: Old repo shows archive banner + notice
3. **Releases**: Old repo shows last release was v6.x, new repo shows v7.0
4. **Stars/Forks**: Users following old repo may see activity in new repo
5. **Documentation**: Your website/docs can point to new version

### Manual Discovery:
- Old repo README prominently displays link to new version
- Archive banner is very visible
- GitHub suggests related repositories

## What Happens to Old Links?

### If Old Repo is Archived (Recommended):
- ✅ `github.com/ranganathanlab/pySCA` → Shows archived repo with notice
- ✅ Clone URLs still work (but point to archived version)
- ✅ Issues/PRs are locked (read-only)
- ✅ Releases are visible but archived
- ✅ Users see clear notice to use new version

### If Old Repo is Private:
- ❌ `github.com/ranganathanlab/pySCA` → 404 or "private repository"
- ❌ Old links break
- ❌ Users can't find new version easily
- ❌ Bad user experience

## Best Practice Recommendation

**Archive, don't privatize:**

1. **Archive old `pySCA`** (v6.x) - keeps it public but read-only
2. **Make `pySCA-dev` public and rename to `pySCA`** (v7.0)
3. **Add clear notices in both repositories**
4. **Update website/documentation** to point to new version

This way:
- ✅ Old links still work
- ✅ Users naturally discover new version
- ✅ Clear migration path
- ✅ Professional approach
- ✅ No broken links or 404s

## Summary

**Question: Will old link point to new release?**
- **If archived**: Old link works, shows archive notice + link to new version
- **If private**: Old link breaks (404), users can't find new version easily

**Question: Will users find it naturally?**
- **If archived**: Yes, via archive notice, search, and repository suggestions
- **If private**: No, they'll get 404 errors

**Recommendation**: Archive the old repository rather than making it private.
