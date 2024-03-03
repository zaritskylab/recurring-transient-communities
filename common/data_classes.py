from dataclasses import dataclass
from typing import Set, List

from mashumaro import DataClassJSONMixin, DataClassDictMixin


class Base(DataClassJSONMixin, DataClassDictMixin):
    pass

@dataclass
class Community(Base):
    cells: Set[int]
    community_id: int
    community_size: int
    total_community_events: int = 1 # complete + partial
    coll_ids: List[int] = None # All coll_ids of this community
    original_coll_ids: List[int] = None # Original coll_id of this community
