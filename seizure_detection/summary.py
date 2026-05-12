import pandas as pd


def build_animal_summary(df_events: pd.DataFrame) -> pd.DataFrame:
    if df_events.empty:
        return pd.DataFrame(
            columns=[
                "group",
                "animal_id",
                "seizure_count",
                "total_duration_sec",
                "mean_duration_sec",
            ]
        )

    summary = (
        df_events.groupby(["group", "animal_id"], as_index=False)
        .agg(
            seizure_count=("event_index", "count"),
            total_duration_sec=("duration_sec", "sum"),
            mean_duration_sec=("duration_sec", "mean"),
        )
    )

    return summary
