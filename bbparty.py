#!/usr/bin/env python

import socket
import errno
from io import BytesIO
from time import sleep, time
import traceback
import ast


class Quit(Exception):
    pass


class Connection(object):

    def get_data(self):
        incoming_size = int(self.socket.recv(16, socket.MSG_WAITALL), 16)
        incoming_data = self.socket.recv(incoming_size, socket.MSG_WAITALL)
        return str(incoming_data, 'utf8')

    def put_data(self, data):
        output = BytesIO()
        output.write(b'\x00' * 16)
        output.write(bytes(data, 'utf8'))
        output_size = b'%016x' % (output.tell() - 16)
        assert len(output_size) == 16
        output.seek(0)
        output.write(output_size)
        output_data = output.getvalue()
        self.socket.sendall(output_data)


class IncomingConnection(Connection):

    def __init__(self, listen_addr, peer_addr):
        self.local_addr = listen_addr
        self.peer_addr = peer_addr
        listen_socket = socket.socket(
            socket.AF_INET,
            socket.SOCK_STREAM,
            socket.IPPROTO_TCP,
        )
        listen_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        listen_socket.settimeout(None)
        listen_socket.bind(self.local_addr)
        listen_socket.listen(1)
        self.listen_socket = listen_socket
        self.socket = None
        self.peer_addr = None
        #print("listening to {0}".format(self.local_addr))

    def initialize(self):
        accepted_socket, accepted_addr = self.listen_socket.accept()
        accepted_socket.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        if self.peer_addr == accepted_addr:
            import pdb; pdb.set_trace()
        self.peer_addr = accepted_addr
        self.socket = accepted_socket
        # print("connection from {0}".format(self.peer_addr))

    def close(self):
        if self.listen_socket is not None:
            try:
                self.listen_socket.shutdown(socket.SHUT_RD)
            except Exception as e:
                traceback.print_exc()
            try:
                self.listen_socket.close()
            except Exception as e:
                traceback.print_exc()

        if self.socket is not None:
            try:
                self.socket.shutdown(socket.SHUT_RD)
            except Exception as e:
                traceback.print_exc()
            try:
                self.socket.close()
            except Exception as e:
                traceback.print_exc()

        self.socket = None


class OutgoingConnection(Connection):

    def __init__(self, connect_addr, local_addr=None):
        self.peer_addr = connect_addr
        self.local_addr = local_addr
        connect_socket = socket.socket(
            socket.AF_INET,
            socket.SOCK_STREAM,
            socket.IPPROTO_TCP,
        )
        connect_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        connect_socket.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        if local_addr is not None:
            connect_socket.bind(local_addr)
        self.socket = connect_socket

    def initialize(self):
        while True:
            try:
                # print("connecting to {0}".format(self.peer_addr))
                self.socket.connect(self.peer_addr)
                break
            except socket.error as e:
                if e.errno == errno.ECONNREFUSED:
                    sleep(1)
                else:
                    raise

    def close(self):
        if self.socket is not None:
            try:
                self.socket.shutdown(socket.SHUT_WR)
            except Exception as e:
                traceback.print_exc()

            try:
                self.socket.close()
            except Exception as e:
                traceback.print_exc()

        self.socket = None


class BBParty(object):

    def __init__(self, peers, name, spec):

        self.peers = peers
        self.base_port = spec.get('base_port', 11000)
        self.spec = spec
        self.connections = []
        self.transcript = []
        self.timestamp = time()
        self.peer_index = self.validate(name)
        self.initialize()

    def validate(self, name):
        peer_index = None
        for i, peer in self.enum_peers():
            if 'host' not in peer:
                m = "peers must all have a 'host' key"
                raise ValueError(m)
            if 'name' not in peer:
                m = "peers must all have a 'name' key"
                raise ValueError(m)
            if name == peer['name']:
                peer_index = i

        if peer_index is None:
            m = "peer {0!r} not found in config"
            m = m.format(name)
            raise KeyError(m)

        return peer_index


    def enum_peers(self):
        return enumerate(self.peers)

    def initialize(self):
        self_host = self.peers[self.peer_index]['host']
        for i, peer in self.enum_peers():
            peer_port = self.base_port + len(self.peers) * i + self.peer_index
            local_port = self.base_port + len(self.peers) * self.peer_index + i
            local_addr = (self_host, local_port)
            peer_addr = (peer['host'], peer_port)
            if i < self.peer_index:
                connection_in = IncomingConnection(local_addr, peer_addr)
                self.connections.append(connection_in)
                connection_in.initialize()
            elif i > self.peer_index:
                connection_out = OutgoingConnection(peer_addr, local_addr)
                self.connections.append(connection_out)
                connection_out.initialize()
            else:
                self.connections.append(None)

    def quit(self):
        for connection in self.connections:
            if connection is None:
                continue
            connection.close()

    def get_step_owner(self, step_index):
        return step_index % len(self.peers)

    def produce_step_output(self):
        # override this
        return len(self.transcript)

    def consume_step_output(self, step_output):
        # override this
        assert len(self.transcript) == step_output
        steps = 10
        if step_output % steps == 0:
            t = time() - self.timestamp
            self.timestamp += t
            m = "output {0!r} sec {1:.3f} steps/sec"
            m = m.format(step_output, steps/t)
            print(m)
        return step_output

    def make_step(self):
        step_index = len(self.transcript)
        step_owner = self.get_step_owner(step_index)
        if step_owner == self.peer_index:
            # make output step
            step_output = self.produce_step_output()
            step_output_data = repr(step_output)
            for connection in self.connections:
                if connection is None:
                    continue
                connection.put_data(step_output_data)
            self.transcript.append(step_output)
        else:
            # make input step
            step_output_data = self.connections[step_owner].get_data()
            step_output = ast.literal_eval(step_output_data)
            accepted_step_output = self.consume_step_output(step_output)
            self.transcript.append(accepted_step_output)

    def run(self):
        try:
            while True:
                self.make_step()
        except Quit as e:
            return


class DemoBBParty(BBParty):
    pass



peers = [
    {
        'name': 'one',
        'host': '127.0.0.1',
        'port': 11001,
    },
    {
        'name': 'two',
        'host': '127.0.0.1',
        'port': 11002,
    },
]


def print_help():
    import sys
    print("Usage: {0} <config_file> <peer_name>".format(sys.argv[0]))
    raise SystemExit(1)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print_help()

    config_file = sys.argv[1]
    with open(config_file) as f:
        config = ast.literal_eval(f.read())

    peers = config['peers']
    spec = config.get('spec', {})

    name = sys.argv[2]
    spec = {}
    bbparty = DemoBBParty(peers, name, spec)
    bbparty.run()
    bbparty.quit()
