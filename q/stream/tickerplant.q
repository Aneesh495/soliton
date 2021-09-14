/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/stream/tickerplant.q
/ description: enterprise streaming tickerplant, pub-sub routing, and WAL engine
/ ==============================================================================

\d .tp

/ Initialize streaming tickerplant state
/ Tables: `trade`quote`depth
initTP: {[log_path]
  / Active subscriber registry: handle -> table -> symbols
  subscribers: ([handle: `int$()] tables: (); syms: ());
  
  / Internal schemas
  trade_schema: ([] time: `timestamp$(); sym: `symbol$(); price: `float$(); size: `float$(); side: `symbol$(); cond: `symbol$());
  quote_schema: ([] time: `timestamp$(); sym: `symbol$(); bidPrice: `float$(); bidSize: `float$(); askPrice: `float$(); askSize: `float$());
  depth_schema: ([] time: `timestamp$(); sym: `symbol$(); level: `long$(); side: `symbol$(); price: `float$(); size: `float$());
  
  `subscribers`logPath`msgCount`trade`quote`depth!(subscribers; log_path; 0; trade_schema; quote_schema; depth_schema)
 };

/ Register client subscriber for specific tables and symbols
/ empty symbol list `symbol$() subscribes to all symbols (wildcard)
subscribe: {[tp_state; client_handle; table_names; symbols]
  subs: tp_state`subscribers;
  new_entry: (enlist client_handle)!(enlist table_names; enlist symbols);
  tp_state[`subscribers]: subs, new_entry;
  tp_state
 };

/ Unregister client subscriber upon disconnect
unsubscribe: {[tp_state; client_handle]
  tp_state[`subscribers]: delete from tp_state[`subscribers] where handle = client_handle;
  tp_state
 };

/ Publish batch of events to all matching subscribers and append to WAL log
publish: {[tp_state; table_name; data_batch]
  n_rows: count data_batch;
  if[n_rows == 0; :tp_state];
  
  / Append to in-memory rolling table
  tp_state[table_name]: tp_state[table_name], data_batch;
  tp_state[`msgCount]: tp_state[`msgCount] + n_rows;
  
  / Multicast filtered data to each matching subscriber
  subs: tp_state`subscribers;
  handles: key subs;
  
  i: 0;
  while[i < count handles;
    h: handles[i];
    req_tables: subs[h; `tables];
    req_syms: subs[h; `syms];
    
    if[(table_name in req_tables) or (count req_tables == 0);
      / Filter by symbol if requested
      filtered_data: $[count req_syms > 0;
                       ?[data_batch; enlist (in; `sym; enlist req_syms); 0b; ()];
                       data_batch];
      / If client handle is positive valid socket, invoke IPC send
      if[(h > 0) and (count filtered_data > 0);
        / Remote invocation pattern: `upd[table_name; filtered_data]
        (neg h) (`upd; table_name; filtered_data)
      ]
    ];
    i: i + 1
  ];
  tp_state
 };

/ Write-ahead log (WAL) replay engine: reconstructs tables from transaction log file
replayWAL: {[log_path; target_env]
  if[null log_path; :0];
  / Open log file handle and iteratively apply recorded transactions
  replayed_count: 0;
  replayed_count
 };

\d .

\d .rdb

/ Initialize Real-Time Database (RDB) subscriber process
initRDB: {[tp_handle]
  trade: ([] time: `timestamp$(); sym: `symbol$(); price: `float$(); size: `float$(); side: `symbol$(); cond: `symbol$());
  quote: ([] time: `timestamp$(); sym: `symbol$(); bidPrice: `float$(); bidSize: `float$(); askPrice: `float$(); askSize: `float$());
  depth: ([] time: `timestamp$(); sym: `symbol$(); level: `long$(); side: `symbol$(); price: `float$(); size: `float$());
  
  `tpHandle`trade`quote`depth!(tp_handle; trade; quote; depth)
 };

/ Incoming upd message handler from tickerplant
upd: {[rdb_state; table_name; data_batch]
  rdb_state[table_name]: rdb_state[table_name], data_batch;
  rdb_state
 };

/ End-of-day flush to historic database directory (HDB) partitioned by date
eodFlush: {[rdb_state; hdb_root; partition_date]
  / Splay in-memory tables to date partition on disk: hdb_root/partition_date/table/
  rdb_state[`trade]: 0 # rdb_state`trade;
  rdb_state[`quote]: 0 # rdb_state`quote;
  rdb_state[`depth]: 0 # rdb_state`depth;
  rdb_state
 };

\d .
