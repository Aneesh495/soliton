#!/usr/bin/env python3
"""
Soliton Synthetic Microstructure Feed Generator
Simulates realistic L2/L3 limit order book events and streaming trades over kdb+ v3 IPC.
"""

import sys
import os
import struct
import socket
import time
import math
import random
import csv

def pack_ipc_header(payload_len, msg_type=0):
    """
    Format 8-byte kdb+ v3 IPC header:
    endian(1) | msg_type(1) | compressed(0) | unused(0) | total_len(4)
    """
    total_len = payload_len + 8
    return struct.pack("<BBBB I", 1, msg_type, 0, 0, total_len)

def generate_synthetic_data(num_events=5000, initial_price=150.0, output_csv=None):
    """
    Generate synthetic order flow governed by a 2D Hawkes self/cross-exciting process.
    Dimension 0: Buy order arrivals
    Dimension 1: Sell order arrivals
    """
    print(f"Generating {num_events} high-frequency market events...")
    
    # Hawkes parameters
    mu = [20.0, 20.0]
    alpha = [[8.0, 4.0],
             [4.0, 8.0]]
    beta = [[15.0, 15.0],
            [15.0, 15.0]]
    
    events = []
    t = 0.0
    mid = initial_price
    spread = 0.05
    
    r_state = [[0.0, 0.0], [0.0, 0.0]]
    
    for i in range(num_events):
        lambda0 = mu[0] + alpha[0][0] * r_state[0][0] + alpha[0][1] * r_state[0][1]
        lambda1 = mu[1] + alpha[1][0] * r_state[1][0] + alpha[1][1] * r_state[1][1]
        total_lambda = max(1.0, lambda0 + lambda1)
        
        dt = random.expovariate(total_lambda)
        t += dt
        
        # Decay recursive state
        for d in range(2):
            for k in range(2):
                r_state[d][k] *= math.exp(-beta[d][k] * dt)
        
        # Decide which channel fired
        prob_buy = lambda0 / total_lambda
        is_buy = random.random() < prob_buy
        fired_type = 0 if is_buy else 1
        
        # Increment jump
        for d in range(2):
            r_state[d][fired_type] += 1.0
            
        # Price drift and micro-movement
        price_shock = (0.01 if is_buy else -0.01) * random.choice([0, 1, 1, 2])
        mid = max(10.0, mid + price_shock)
        
        best_bid = round(mid - 0.5 * spread, 2)
        best_ask = round(mid + 0.5 * spread, 2)
        
        trade_size = random.choice([100, 200, 300, 500, 1000])
        trade_price = best_ask if is_buy else best_bid
        side_str = "buy" if is_buy else "sell"
        
        events.append({
            "time_nanos": int(t * 1e9),
            "sym": "AAPL",
            "price": trade_price,
            "size": float(trade_size),
            "side": side_str,
            "bestBid": best_bid,
            "bestAsk": best_ask,
            "midPrice": mid
        })
    
    if output_csv:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(events[0].keys()))
            writer.writeheader()
            writer.writerows(events)
        print(f"Successfully saved {num_events} records to {output_csv}")
        
    return events

def stream_to_ipc_server(events, host="127.0.0.1", port=5001):
    """
    Connect to Soliton IPC server and stream events over raw TCP socket.
    """
    print(f"Connecting to Soliton IPC server at {host}:{port}...")
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
        print("Connected! Streaming market events...")
        for ev in events:
            # Build q command expression
            expr = f"trade: trade, enlist `time`sym`price`size`side!({ev['time_nanos']}; `AAPL; {ev['price']}; {ev['size']}; `{ev['side']})"
            payload = expr.encode("utf-8")
            # Wrap as string atom / IPC byte message
            hdr = pack_ipc_header(len(payload) + 1, msg_type=0)
            msg = hdr + b"\x0a" + payload # 10 = char vector
            s.sendall(msg)
        s.close()
        print("Stream completed successfully.")
    except ConnectionRefusedError:
        print(f"Could not connect to {host}:{port}. Run './bin/soliton -p {port}' first.")

def main():
    csv_path = os.path.join(os.path.dirname(__file__), "..", "data", "sample_ticks.csv")
    events = generate_synthetic_data(num_events=5000, output_csv=csv_path)
    if "--stream" in sys.argv:
        stream_to_ipc_server(events)

if __name__ == "__main__":
    main()
