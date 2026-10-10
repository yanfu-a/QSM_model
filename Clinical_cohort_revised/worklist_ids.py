"""Consistent handling of repeated image_dir_id values in cohort worklists."""


def resolve_duplicate_ids(worklist):
    """Keep identical duplicate rows once; exclude every ID with conflicting rows.

    Check the complete worklist before applying runnable, model or ID filters. Otherwise
    one of two conflicting rows may be filtered away and the remaining row accepted.
    Returns (unique_rows, exclusions), with one exclusion record per conflicting ID.
    """
    repeated = worklist[worklist["image_dir_id"].duplicated(keep=False)]
    conflicts = []
    for oid, group in repeated.groupby("image_dir_id", sort=False, dropna=False):
        differing = [
            f"{column}={' | '.join(group[column].astype(str).unique())}"
            for column in worklist.columns if column != "image_dir_id"
            and group[column].nunique(dropna=False) > 1
        ]
        if differing:
            conflicts.append({"IID": oid, "reason": "worklist rows disagree: " + "; ".join(differing)})
    rejected = {record["IID"] for record in conflicts}
    unique_rows = worklist[~worklist["image_dir_id"].isin(rejected)].drop_duplicates("image_dir_id")
    return unique_rows, conflicts
